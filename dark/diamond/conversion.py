import bz2
from json import dumps, loads
from operator import itemgetter
from collections import Counter

from Bio.File import as_handle

from dark.hsp import HSP, LSP
from dark.score import HigherIsBetterScore
from dark.alignments import Alignment, ReadAlignments
from dark.diamond.hsp import normalizeHSP

# The following are the fields (in the order they are expected on the
# command line via --outfmt 6) that must be given to DIAMOND blastx to
# allow its output to be parsed by convert-diamond-to-sam.py (which uses
# diamondTabularFormatToDicts (below)).
FIELDS = (
    "bitscore btop qframe qend full_qqual qlen full_qseq qseqid "
    "qstart slen sstart stitle"
)

# The keys in the following are DIAMOND format 6 field names. The values
# are one-argument functions that take a string and return an appropriately
# converted field value.
#
# The following fields are taken from the DIAMOND manual v0.9.22 2018-05-11.
# Fields whose name doesn't appear here will be left as strings.
DIAMOND_FIELD_CONVERTER = {
    "bitscore": float,
    "evalue": float,
    "frame": int,
    "gapopen": int,
    "gaps": int,
    "identicalCount": lambda nident: None if nident is None else int(nident),
    "length": int,
    "mismatch": int,
    "nident": int,
    "pident": float,
    "positive": int,
    "positiveCount": lambda pos: None if pos is None else int(pos),
    "ppos": float,
    "qcovhsp": float,
    "qend": int,
    "qframe": int,
    "qlen": int,
    "qstart": int,
    "score": float,
    "send": int,
    "slen": int,
    "sstart": int,
}


class DiamondTabularFormat:
    """
    Read/convert DIAMOND TAB-separated format (format #6).

    @param fieldNames: A C{list} or C{tuple} of C{str} DIAMOND field names.
        Run 'diamond -help' to see the full list. If C{None}, a default set of
        fields will be used, as compatible with convert-diamond-to-sam.py
    """

    def __init__(self, fieldNames=None):
        self._fieldNames = fieldNames or FIELDS.split()
        self._nFields = len(self._fieldNames)
        if not self._nFields:
            raise ValueError("fieldNames cannot be empty.")

        c = Counter(self._fieldNames)
        if c.most_common(1)[0][1] > 1:
            raise ValueError(
                "field names contains duplicated names: %s."
                % (", ".join(sorted(x[0] for x in c.most_common() if x[1] > 1)))
            )

        def identity(x):
            return x

        self._convertors = tuple(
            DIAMOND_FIELD_CONVERTER.get(field, identity) for field in self._fieldNames
        )

    def diamondFieldsToDict(self, line):
        """
        Convert a line of DIAMOND tabular (--outfmt 6) output to a dictionary

        @param line: A C{str} line read from DIAMOND output. This may contain
            a trailing whitespace (including newline and carriage return).
        @raise ValueError: If a line of C{filename} does not have the expected
            number of TAB-separated fields (i.e., len(self._fieldNames)).
        @return: A C{dict} with keys that are the DIAMOND field names, with
            values as converted by the functions set up in __init__.
        """
        values = line.rstrip().split("\t")

        if len(values) != self._nFields:
            raise ValueError(
                "DIAMOND output line had %d field values (expected %d). "
                "The offending input line was %r." % (len(values), self._nFields, line)
            )

        return dict(
            (fieldName, func(value))
            for fieldName, func, value in zip(
                self._fieldNames, self._convertors, values
            )
        )

    def diamondTabularFormatToDicts(self, filename):
        """
        Read DIAMOND tabular (--outfmt 6) output and convert lines to
        dictionaries.

        @param filename: Either a C{str} file name or an open file pointer.
        @raise ValueError: If a line of C{filename} does not have the expected
            number of TAB-separated fields (i.e., len(fieldNames)).
        @return: A generator that yields C{dict}s returned by
            C{diamondFieldsToDict}, above.
        """
        with as_handle(filename) as fp:
            for line in fp:
                yield self.diamondFieldsToDict(line)


class DiamondTabularFormatReader:
    """
    Provide a method that yields parsed tabular records from a file. Store and
    make accessible the global DIAMOND parameters.

    Make sure you run DIAMOND with the right output format. You must have run
    with --outfmt 6 and one of the following sets of arguments (in the given
    order):

        13 arg version:
        qtitle stitle bitscore evalue qframe qseq qstart qend sseq sstart send
            slen btop

        15 arg version:
        qtitle stitle bitscore evalue qframe qseq qstart qend sseq sstart send
            slen btop nident positive

        17 arg version:
        qtitle stitle bitscore evalue qframe qseq qstart qend sseq sstart send
            slen btop nident pident positive ppos

    @param filename: A C{str} filename or an open file pointer, containing
        DIAMOND tabular records.

    """

    def __init__(self, filename):
        self._filename = filename
        self.application = "DIAMOND"
        self.params = {
            "application": self.application,
            "reference": (
                "Buchfink et al., Fast and Sensitive Protein "
                "Alignment using DIAMOND, Nature Methods, 12, 59-60 "
                "(2015)"
            ),
            "task": "blastx",  # TODO: Add support for blastp, if needed.
            "version": "v0.8.23",
        }

    def records(self):
        """
        Parse the DIAMOND output and yield records. This will be used to read
        original DIAMOND output (either from stdin or from a file) to turn the
        DIAMOND results into Python dictionaries that will then be stored in
        our JSON format.

        @return: A generator that produces C{dict}s containing 'alignments' and
            'query' C{str} keys.
        """
        with as_handle(self._filename) as fp:
            previousQtitle = None
            subjectsSeen = {}
            record = {}
            for line in fp:
                line = line[:-1]
                fields = line.split("\t")
                nFields = len(fields)

                if nFields == 13:
                    (
                        qtitle,
                        stitle,
                        bitscore,
                        evalue,
                        qframe,
                        qseq,
                        qstart,
                        qend,
                        sseq,
                        sstart,
                        send,
                        slen,
                        btop,
                    ) = fields
                    nident = pident = positive = ppos = None
                elif nFields == 15:
                    (
                        qtitle,
                        stitle,
                        bitscore,
                        evalue,
                        qframe,
                        qseq,
                        qstart,
                        qend,
                        sseq,
                        sstart,
                        send,
                        slen,
                        btop,
                        nident,
                        positive,
                    ) = fields
                    pident = ppos = None
                elif nFields == 17:
                    (
                        qtitle,
                        stitle,
                        bitscore,
                        evalue,
                        qframe,
                        qseq,
                        qstart,
                        qend,
                        sseq,
                        sstart,
                        send,
                        slen,
                        btop,
                        nident,
                        pident,
                        positive,
                        ppos,
                    ) = fields
                else:
                    raise ValueError(
                        "Could not make sense of DIAMOND output. You must use "
                        "--outfmt 6 with 13, 15, or 17 arguments. See %s for "
                        "the expected names and order." % __file__
                    )

                hsp = {
                    "bits": float(bitscore),
                    "btop": btop,
                    "expect": float(evalue),
                    "frame": int(qframe),
                    "identicalCount": None if nident is None else int(nident),
                    "percentIdentical": (None if pident is None else float(pident)),
                    "positiveCount": (None if positive is None else int(positive)),
                    "percentPositive": (None if ppos is None else float(ppos)),
                    "query": qseq,
                    "query_start": int(qstart),
                    "query_end": int(qend),
                    "sbjct": sseq,
                    "sbjct_start": int(sstart),
                    "sbjct_end": int(send),
                }

                if previousQtitle == qtitle:
                    # We have already started accumulating alignments for this
                    # query.
                    if stitle not in subjectsSeen:
                        # We have not seen this subject before, so this is a
                        # new alignment.
                        subjectsSeen.add(stitle)
                        alignment = {
                            "hsps": [hsp],
                            "length": int(slen),
                            "title": stitle,
                        }
                        record["alignments"].append(alignment)
                    else:
                        # We have already seen this subject, so this is another
                        # HSP in an already existing alignment.
                        for alignment in record["alignments"]:
                            if alignment["title"] == stitle:
                                alignment["hsps"].append(hsp)
                                break
                else:
                    # All alignments for the previous query id (if any)
                    # have been seen.
                    if previousQtitle is not None:
                        yield record

                    # Start building up the new record.
                    record = {}
                    subjectsSeen = {stitle}
                    alignment = {
                        "hsps": [hsp],
                        "length": int(slen),
                        "title": stitle,
                    }
                    record["alignments"] = [alignment]
                    record["query"] = qtitle

                    previousQtitle = qtitle

            # Yield the last record, if any.
            if record:
                yield record

    def saveAsJSON(self, fp, writeBytes=False):
        """
        Write the records out as JSON. The first JSON object saved contains
        information about the DIAMOND algorithm.

        @param fp: A C{str} file pointer to write to.
        @param writeBytes: If C{True}, the JSON will be written out as bytes
            (not strings). This is required when we are writing to a BZ2 file.
        """
        if writeBytes:
            fp.write(dumps(self.params, sort_keys=True).encode("UTF-8"))
            fp.write(b"\n")
            for record in self.records():
                fp.write(dumps(record, sort_keys=True).encode("UTF-8"))
                fp.write(b"\n")
        else:
            fp.write(dumps(self.params, sort_keys=True))
            fp.write("\n")
            for record in self.records():
                fp.write(dumps(record, sort_keys=True))
                fp.write("\n")


class JSONRecordsReader:
    """
    Provide a method that yields JSON records from a file. Store, check, and
    make accessible the DIAMOND parameters.

    @param filename: A C{str} filename containing JSON DIAMOND records.
    @param scoreClass: A class to hold and compare scores (see scores.py).
        Default is C{HigherIsBetterScore}, for comparing bit scores. If you
        are using e-values, pass LowerIsBetterScore instead.
    """

    def __init__(self, filename, scoreClass=HigherIsBetterScore):
        self._filename = filename
        self._scoreClass = scoreClass
        if scoreClass is HigherIsBetterScore:
            self._hspClass = HSP
        else:
            self._hspClass = LSP

        self._open(filename)
        self.diamondTask = self.params["task"]

    def _open(self, filename):
        """
        Open the input file. Set self._fp to point to it. Read the first
        line of parameters.

        @param filename: A C{str} filename containing JSON DIAMOND records.
        @raise ValueError: if the first line of the file isn't valid JSON,
            if the input file is empty, or if the JSON does not contain an
            'application' key.
        """
        if filename.endswith(".bz2"):
            self._fp = bz2.open(filename, mode="rt", encoding="UTF-8")
        else:
            self._fp = open(filename)

        line = self._fp.readline()
        if not line:
            raise ValueError("JSON file %r was empty." % self._filename)

        try:
            self.params = loads(line[:-1])
        except ValueError as e:
            raise ValueError(
                "Could not convert first line of %r to JSON (%s). "
                "Line is %r." % (self._filename, e, line[:-1])
            )

    def _dictToAlignments(self, diamondDict, read):
        """
        Take a dict (made by DiamondTabularFormatReader.records)
        and convert it to a list of alignments.

        @param diamondDict: A C{dict}, from records().
        @param read: A C{Read} instance, containing the read that DIAMOND used
            to create this record.
        @return: A C{list} of L{dark.alignment.Alignment} instances.
        """
        alignments = []
        getScore = itemgetter("bits" if self._hspClass is HSP else "expect")

        for diamondAlignment in diamondDict["alignments"]:
            alignment = Alignment(diamondAlignment["length"], diamondAlignment["title"])
            alignments.append(alignment)
            for diamondHsp in diamondAlignment["hsps"]:
                score = getScore(diamondHsp)
                normalized = normalizeHSP(diamondHsp, len(read), self.diamondTask)
                hsp = self._hspClass(
                    score,
                    readStart=normalized["readStart"],
                    readEnd=normalized["readEnd"],
                    readStartInSubject=normalized["readStartInSubject"],
                    readEndInSubject=normalized["readEndInSubject"],
                    readFrame=diamondHsp["frame"],
                    subjectStart=normalized["subjectStart"],
                    subjectEnd=normalized["subjectEnd"],
                    readMatchedSequence=diamondHsp["query"],
                    subjectMatchedSequence=diamondHsp["sbjct"],
                    # Use diamondHsp.get on identicalCount, positiveCount,
                    # percentPositive, and percentIdentical because they
                    # were either added in version 2.0.3 or we didn't start
                    # using them until much later and so will not be
                    # present in any of our JSON output generated before
                    # that. Those values will be None when reading those
                    # JSON files, but that's much better than no longer
                    # being able to read all that earlier data.
                    identicalCount=diamondHsp.get("identicalCount"),
                    positiveCount=diamondHsp.get("positiveCount"),
                    percentIdentical=diamondHsp.get("percentIdentical"),
                    percentPositive=diamondHsp.get("percentPositive"),
                )

                alignment.addHsp(hsp)

        return alignments

    def readAlignments(self, reads):
        """
        Read lines of JSON from self._filename, convert them to read alignments
        and yield them.

        @param reads: An iterable of L{Read} instances, corresponding to the
            reads that were given to DIAMOND.
        @raise ValueError: If any of the lines in the file cannot be converted
            to JSON.
        @return: A generator that yields C{dark.alignments.ReadAlignments}
            instances.
        """
        if self._fp is None:
            self._open(self._filename)

        reads = iter(reads)

        try:
            for lineNumber, line in enumerate(self._fp, start=2):
                try:
                    record = loads(line[:-1])
                except ValueError as e:
                    raise ValueError(
                        "Could not convert line %d of %r to JSON (%s). "
                        "Line is %r." % (lineNumber, self._filename, e, line[:-1])
                    )
                else:
                    recordTitle = record["query"]
                    while True:
                        # Iterate through the input reads until we find the
                        # one that matches this DIAMOND record.
                        try:
                            read = next(reads)
                        except StopIteration:
                            raise ValueError(
                                "Read generator failed to yield a read "
                                "with id '%s' as found in record number %d "
                                "during parsing of DIAMOND output file %r."
                                % (recordTitle, lineNumber - 1, self._filename)
                            )
                        else:
                            # Look for an exact read id / subject title match.
                            # If that doesn't work, allow for the case where
                            # the JSON record has a truncated query (i.e.,
                            # read) id. This covers the situation where a tool
                            # we use (e.g., bwa mem) unconditionally does this
                            # truncation in the output it writes.
                            if (
                                read.id == recordTitle
                                or read.id.split()[0] == recordTitle
                            ):
                                alignments = self._dictToAlignments(record, read)
                                yield ReadAlignments(read, alignments)
                                break
                            else:
                                # This is an input read that had no DIAMOND
                                # matches. So it does not appear in the
                                # DIAMOND's output. Yield an empty
                                # ReadAlignments for it.
                                yield ReadAlignments(read, [])

        finally:
            self._fp.close()
            self._fp = None
