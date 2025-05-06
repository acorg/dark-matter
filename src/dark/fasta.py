from hashlib import md5
import sqlite3
import os
from typing import Generator, Iterator, Iterable, TextIO, Union, Type, Optional, List

from Bio import SeqIO, bgzf  # type: ignore

from dark.reads import Reads, DNARead, Read
from dark.sqlite3 import sqliteConnect
from dark.utils import asHandle


def fastaToList(fastaFilename) -> list:
    return list(SeqIO.parse(fastaFilename, "fasta"))


def dedupFasta(reads: Reads) -> Generator[Read, None, None]:
    """
    Remove sequence duplicates (based on sequence) from FASTA.

    @param reads: a C{dark.reads.Reads} instance.
    @return: a generator of C{dark.reads.Read} instances with no duplicates.
    """
    seen: set[bytes] = set()
    add = seen.add

    for read in iter(reads):
        hash_ = md5(read.sequence.encode("UTF-8")).digest()
        if hash_ not in seen:
            add(hash_)
            yield read


def dePrefixAndSuffixFasta(
    sequences: Iterable[SeqIO.SeqRecord],
) -> Iterator[SeqIO.SeqRecord]:
    """
    sequences: an iterator producing Bio.Seq sequences.

    @return: a generator of sequences with no duplicates and no fully contained
        subsequences.
    """
    sequences = sorted(sequences, key=lambda s: len(s.seq), reverse=True)
    seen = set()
    for s in sequences:
        thisSeq = str(s.seq)
        thisHash = md5(thisSeq.encode("UTF-8")).digest()
        if thisHash not in seen:
            # Add prefixes.
            newHash = md5()
            for nucl in thisSeq:
                newHash.update(nucl.encode("UTF-8"))
                seen.add(newHash.digest())
            # Add suffixes.
            for start in range(len(thisSeq) - 1):
                seen.add(md5(thisSeq[start + 1 :].encode("UTF-8")).digest())
            yield s


def fastaSubtract(fastaFiles: List[TextIO]) -> Iterator[SeqIO.SeqRecord]:
    """
    Given a list of open file descriptors, each with FASTA content,
    remove the reads found in the 2nd, 3rd, etc files from the first file
    in the list.

    @param fastaFiles: a C{list} of FASTA filenames.
    @raises IndexError: if passed an empty list.
    @return: An iterator producing C{Bio.Seq} instances suitable for
        writing to a file using C{Bio.SeqIO.write}.

    """
    reads = {}
    fastaFiles = list(fastaFiles)
    firstFile = fastaFiles.pop(0)
    for seq in SeqIO.parse(firstFile, "fasta"):
        reads[seq.id] = seq

    for fastaFile in fastaFiles:
        for seq in SeqIO.parse(fastaFile, "fasta"):
            # Make sure that reads with the same id have the same sequence.
            if seq.id in reads:
                assert str(seq.seq) == str(reads[seq.id].seq)
            reads.pop(seq.id, None)

    return iter(reads.values())


class FastaReads(Reads):
    """
    Subclass of L{dark.reads.Reads} providing access to FASTA reads.

    @param _files: Either a single C{str} file name or file handle, or a
        C{list} of C{str} file names and/or file handles. Each file or file
        handle must contain sequences in FASTA format.
    @param readClass: The class of read that should be yielded by iter.
    @param upperCase: If C{True}, read sequences will be converted to upper
        case.
    """

    def __init__(
        self,
        _files: Union[str, TextIO, Union[list, tuple]],
        readClass: Type[Read] = DNARead,
        upperCase: bool = False,
    ):
        self._files = _files if isinstance(_files, (list, tuple)) else [_files]
        self._readClass = readClass
        # TODO: It would be better if upperCase were an argument that could
        # be passed to Reads.__init__ and that could do the uppercasing in
        # its add method (as opposed to using it below in our iter method).
        # In that case, in the iter of this class we'd call self.add on
        # each of the sequences coming from self._file. Or, if we'd already
        # read the file we'd return Reads.iter(self) to re-iterate over the
        # sequences already added from the file.
        self._upperCase = upperCase
        super().__init__()

    def iter(self) -> Generator[Read, None, None]:
        """
        Iterate over the sequences in the files in self.files_, yielding each
        as an instance of the desired read class.
        """
        for _file in self._files:
            with asHandle(_file) as fp:
                # Duplicate some code here so as not to test
                # self._upperCase in the loop.
                if self._upperCase:
                    for seq in SeqIO.parse(fp, "fasta"):
                        read = self._readClass(seq.description, str(seq.seq.upper()))
                        yield read
                else:
                    for seq in SeqIO.parse(fp, "fasta"):
                        read = self._readClass(seq.description, str(seq.seq))
                        yield read


def combineReads(
    filename: str,
    sequences: list[str],
    readClass: Type[Read] = DNARead,
    upperCase: bool = False,
    idPrefix: str = "command-line-read-",
) -> Union[FastaReads, Reads]:
    """
    Combine FASTA reads from a file and/or sequence strings.

    @param filename: A C{str} file name containing FASTA reads.
    @param sequences: A C{list} of C{str} sequences. If a sequence
        contains spaces, the last field (after splitting on spaces) will be
        used as the sequence and the first fields will be used as the sequence
        id.
    @param readClass: The class of the individual reads.
    @param upperCase: If C{True}, reads will be converted to upper case.
    @param idPrefix: The C{str} prefix that will be used for the id of the
        sequences in C{sequences} that do not have an id specified. A trailing
        sequence number will be appended to this prefix. Note that
        'command-line-read-', the default id prefix, could collide with ids in
        the FASTA file, if given. So output might be ambiguous. That's why we
        allow the caller to specify a custom prefix.
    @return: A C{FastaReads} instance.
    """
    # Read sequences from a FASTA file, if given.
    reads: Union[FastaReads, Reads]

    if filename:
        reads = FastaReads(filename, readClass=readClass, upperCase=upperCase)
    else:
        reads = Reads()

    # Add any individually specified subject sequences.
    if sequences:
        for count, sequence in enumerate(sequences, start=1):
            # Try splitting the sequence on its last space and using the
            # first part of the split as the read id. If there's no space,
            # assign a generic id.
            parts = sequence.rsplit(" ", 1)
            if len(parts) == 2:
                readId, sequence = parts
            else:
                readId = "%s%d" % (idPrefix, count)
            if upperCase:
                sequence = sequence.upper()
            read = readClass(readId, sequence)
            reads.add(read)

    return reads


class SqliteIndex:
    """
    Create an Sqlite3 database holding FASTA sequence ids, file names, and
    offsets for fast random dictionary-like access.

    @param dbFilename: A C{str} file name containing an sqlite3 database. If
        the file does not exist it will be created. The special string
        ":memory:" can be used to create an in-memory database.
    @param readClass: The class of read that should be returned by __getitem__.
    @param fastaDirectory: A C{str} directory where the indexed FASTA files
        can be found. If provided, this directory is only used by __getitem__,
        which will combine it with the basename of the files given to
        C{addFile} to locate the FASTA.
    """

    def __init__(
        self,
        dbFilename: str,
        readClass: Type[Read] = DNARead,
        fastaDirectory: Optional[str] = None,
    ):
        self._readClass = readClass
        self._fastaDirectory = fastaDirectory
        creating = dbFilename == ":memory:" or not os.path.exists(dbFilename)
        self._connection: sqlite3.Connection = sqliteConnect(dbFilename)
        if creating:
            # Create a new database.
            cur = self._connection.cursor()
            cur.executescript(
                """
                CREATE TABLE files (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name VARCHAR UNIQUE
                );

                CREATE TABLE sequences (
                    id VARCHAR UNIQUE PRIMARY KEY,
                    fileNumber INTEGER,
                    offset INTEGER
                );
            """
            )
            self._connection.commit()

    def _getFilename(self, fileNumber: int) -> Optional[str]:
        """
        Given a file number, get its name (if any).

        @param fileNumber: An C{int} file number.
        @return: A C{str} file name or C{None} if a file with that number
            has not been added.
        """
        cur = self._connection.cursor()
        cur.execute("SELECT name FROM files WHERE id = ?", (fileNumber,))
        row = cur.fetchone()
        if row is None:
            return None
        else:
            return row[0]

    def _getFileNumber(self, filename: str) -> Optional[int]:
        """
        Given a file name, get its file number (if any).

        @param filename: A C{str} file name.
        @return: An C{int} file number or C{None} if no file with that name
            has been added.
        """
        cur = self._connection.cursor()
        cur.execute("SELECT id FROM files WHERE name = ?", (filename,))
        row = cur.fetchone()
        if row is None:
            return None
        else:
            return row[0]

    def _addFilename(self, filename: str) -> Optional[int]:
        """
        Add a new file name.

        @param filename: A C{str} file name.
        @raise ValueError: If a file with this name has already been added.
        @return: The C{int} id of the newly added file.
        """
        cur = self._connection.cursor()
        try:
            cur.execute("INSERT INTO files(name) VALUES (?)", (filename,))
        except sqlite3.IntegrityError as e:
            if str(e).find("UNIQUE constraint failed") > -1:
                raise ValueError("Duplicate file name: %r" % filename)
            else:
                raise
        else:
            fileNumber = cur.lastrowid
            self._connection.commit()
            return fileNumber

    def addFile(self, filename: str) -> int:
        """
        Add a new FASTA file of sequences.

        @param filename: A C{str} file name, with the file in FASTA format.
            This file must (obviously) exist at indexing time. When __getitem__
            is used to access sequences, it is possible to provide a
            C{fastaDirectory} argument to our C{__init__} to indicate the
            directory containing the original FASTA files, in which case the
            basename of the file here provided in C{filename} is used to find
            the file in the given directory. This allows the construction of a
            sqlite database from the shell in one directory and its use
            programmatically from another directory.
        @raise ValueError: If a file with this name has already been added or
            if the file contains a sequence whose id has already been seen.
        @return: The C{int} number of sequences added from the file.
        """
        endswith = filename.lower().endswith
        if endswith(".bgz") or endswith(".gz"):
            useBgzf = True
        elif endswith(".bz2"):
            raise ValueError(
                "Compressed FASTA is only supported in BGZF format. Use "
                "bgzip to compresss your FASTA."
            )
        else:
            useBgzf = False

        fileNumber = self._addFilename(filename)
        connection = self._connection
        count = 0
        try:
            with connection:
                if useBgzf:
                    try:
                        fp = bgzf.open(filename, "rb")
                    except ValueError as e:
                        if str(e).find("BGZF") > -1:
                            raise ValueError(
                                "Compressed FASTA is only supported in BGZF "
                                "format. Use the samtools bgzip utility "
                                "(instead of gzip) to compresss your FASTA."
                            )
                        else:
                            raise
                    else:
                        try:
                            for line in fp:
                                if line[0] == ">":
                                    count += 1
                                    id_ = line[1:].rstrip(" \t\n\r")
                                    connection.execute(
                                        "INSERT INTO sequences(id, "
                                        "fileNumber, offset) VALUES (?, ?, ?)",
                                        (id_, fileNumber, fp.tell()),
                                    )
                        finally:
                            fp.close()
                else:
                    with open(filename) as fp:
                        offset = 0
                        for line in fp:
                            offset += len(line)
                            if line[0] == ">":
                                count += 1
                                id_ = line[1:].rstrip(" \t\n\r")
                                connection.execute(
                                    "INSERT INTO sequences(id, fileNumber, "
                                    "offset) VALUES (?, ?, ?)",
                                    (id_, fileNumber, offset),
                                )
        except sqlite3.IntegrityError as e:
            if str(e).find("UNIQUE constraint failed") > -1:
                original = self._find(id_)
                if original is None:
                    # The id must have appeared twice in the current file,
                    # because we could not look it up in the database
                    # (i.e., it was INSERTed but not committed).
                    raise ValueError(
                        "FASTA sequence id '%s' found twice in file '%s'."
                        % (id_, filename)
                    )
                else:
                    origFilename, _ = original
                    raise ValueError(
                        "FASTA sequence id '%s', found in file '%s', was "
                        "previously added from file '%s'."
                        % (id_, filename, origFilename)
                    )
            else:
                raise
        else:
            return count

    def _find(self, id_: str) -> Optional[tuple[Optional[str], int]]:
        """
        Find the filename and offset of a sequence, given its id.

        @param id_: A C{str} sequence id.
        @return: A 2-tuple, containing the C{str} file name and C{int} offset
            within that file of the sequence.
        """
        cur = self._connection.cursor()
        cur.execute("SELECT fileNumber, offset FROM sequences WHERE id = ?", (id_,))
        row = cur.fetchone()
        if row is None:
            return None
        else:
            return self._getFilename(row[0]), row[1]

    def __getitem__(self, id_: str) -> Read:
        """
        Return a read, given its id.

        @param id_: A C{str} sequence id.
        @raise KeyError: If C{id_} is not a known sequence.
        @return: A read of our read class.
        """
        location = self._find(id_)
        if location is None:
            raise KeyError("Unknown sequence: %r" % id_)
        else:
            filename, offset = location
            # If a FASTA directory was provided, look for the FASTA files
            # there, otherwise use the filename that was given to addFile.
            assert filename
            if self._fastaDirectory:
                filename = os.path.join(
                    self._fastaDirectory, os.path.basename(filename)
                )

            endswith = filename.lower().endswith
            if endswith(".bgz") or endswith(".gz"):
                opener = bgzf.open
            else:
                opener = open

            sequence = ""
            with opener(filename) as fp:
                fp.seek(offset)
                while True:
                    line = fp.readline()
                    if not line:
                        # EOF
                        break
                    elif line[0] == ">":
                        # We found the next sequence identifier.
                        break
                    else:
                        sequence += line.rstrip("\n\r")

            return self._readClass(id_, sequence)

    def close(self) -> None:
        self._connection.close()
