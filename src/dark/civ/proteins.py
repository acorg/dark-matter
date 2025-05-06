import os
import re
import sqlite3
import sys
import numpy as np
from Bio import SeqIO
from cachetools import LRUCache, cachedmethod
from collections import defaultdict
from functools import partial
from json import load
from operator import attrgetter, itemgetter
from os.path import dirname, exists, join
from urllib.parse import quote
from textwrap import fill

from dark.dimension import dimensionalIterator
from dark.errors import DatabaseDuplicationError
from dark.fasta import FastaReads
from dark.fastq import FastqReads
from dark.filter import TitleFilter
from dark.genbank import getCDSInfo, getSourceInfo
from dark.html import NCBISequenceLinkURL, NCBISequenceLink, readCountText
from dark.reads import Reads
from dark.sqlite3 import sqliteConnect
from dark.taxonomy import (
    lineageTaxonomyLinks,
    Hierarchy,
    LineageElement,
    isAllowedTaxonomicRank,
)
from dark.utils import asHandle


class PathogenSampleFiles:
    """
    Maintain a cache of FASTA/FASTQ file names for the samples that contain a
    given pathogen, create de-duplicated (by read id) FASTA/FASTQ files
    for each pathogen/sample pair, provide functions to write out index files
    of samples numbers (which are generated here in C{self.add}),
    and provide a filename lookup function for pathogen/sample combinations
    or just pathogen accessions by themselves.

    @param proteinGrouper: An instance of C{ProteinGrouper}.
    @param format_: A C{str}, either 'fasta' or 'fastq' indicating the format
        of the files containing the reads matching proteins.
    @raise ValueError: If C{format_} is unknown.
    """

    def __init__(self, proteinGrouper, format_="fasta"):
        self._proteinGrouper = proteinGrouper
        if format_ in ("fasta", "fastq"):
            self._format = format_
            self._readsClass = FastaReads if format_ == "fasta" else FastqReads
        else:
            raise ValueError("format_ must be either 'fasta' or 'fastq'.")
        self._pathogens = {}
        self._samples = {}
        self._readsFilenames = {}

    def add(self, genomeAccession, sampleName):
        """
        Add a (pathogen accession number, sample name) combination and get its
        FASTA/FASTQ file name and unique read count. Write the FASTA/FASTQ file
        if it does not already exist. Save the unique read count into
        C{self._proteinGrouper}.

        @param genomeAccession: A C{str} pathogen accession number.
        @param sampleName: A C{str} sample name.
        @return: A C{str} giving the FASTA/FASTQ file name holding all the
            reads (without duplicates, by id) from the sample that matched the
            proteins in the given pathogen.
        """
        sampleIndex = self._samples.setdefault(sampleName, len(self._samples))

        try:
            return self._readsFilenames[(genomeAccession, sampleIndex)]
        except KeyError:
            reads = Reads()
            for proteinMatch in self._proteinGrouper.genomeAccessions[genomeAccession][
                sampleName
            ]["proteins"].values():
                for read in self._readsClass(proteinMatch["readsFilename"]):
                    reads.add(read)
            saveFilename = join(
                proteinMatch["outDir"],
                "pathogen-%s-sample-%d.%s"
                % (genomeAccession, sampleIndex, self._format),
            )
            reads.filter(removeDuplicatesById=True)
            nReads = reads.save(saveFilename, format_=self._format)
            # Save the unique read count into self._proteinGrouper
            self._proteinGrouper.genomeAccessions[genomeAccession][sampleName][
                "uniqueReadCount"
            ] = nReads
            self._readsFilenames[(genomeAccession, sampleIndex)] = saveFilename
            return saveFilename

    def lookup(self, genomeAccession, sampleName):
        """
        Look up a pathogen accession number, sample name combination and get
        its FASTA/FASTQ file name.

        This method should be used instead of C{add} in situations where
        you want an exception to be raised if a pathogen/sample combination has
        not already been passed to C{add}.

        @param genomeAccession: A C{str} pathogen accession number.
        @param sampleName: A C{str} sample name.
        @raise KeyError: If the pathogen accession number or sample name have
            not been seen, either individually or in combination.
        @return: A C{str} filename retrieved from self._readsFilenames
        """
        return self._readsFilenames[(genomeAccession, self._samples[sampleName])]

    def writeSampleIndex(self, fp):
        """
        Write a file of sample indices and names, sorted by index.

        @param fp: A file-like object, opened for writing.
        """
        print(
            "\n".join(
                "%d %s" % (index, name)
                for (index, name) in sorted(
                    (index, name) for (name, index) in self._samples.items()
                )
            ),
            file=fp,
        )


class ProteinGrouper:
    """
    Group matched proteins by the pathogen they come from.

    @param proteinGenomeDatabase: A connection to an Sqlite3 database
        holding protein and genome information, as built by
        C{make-protein-database.py}.
    @param assetDir: The C{str} directory name where
        C{noninteractive-alignment-panel.py} put its HTML, blue plot and
        alignment panel images, and FASTA or FASTQ files. This must be relative
        to the filenames that will later be passed to C{addFile}.
    @param sampleName: A C{str} sample name. This takes precedence over
        C{sampleNameRegex} (the two cannot be used together, obviously).
    @param sampleNameRegex: A C{str} regular expression that can be used to
        extract a short sample name from full file names subsequently passed
        to C{self.addFile}. The regular expression must have a matching group
        (delimited by parentheses) to capture the part of the file name that
        should be used as the sample name.
    @param format_: A C{str}, either 'fasta' or 'fastq' indicating the format
        of the files containing the reads matching proteins.
    @param saveReadLengths: If C{True}, save the lengths of all reads matching
        proteins.
    @param titleRegex: A regex that pathogen names must match.
        Note that this matching is done on the final part of the protein title
        in square brackets, according to the convention used by the NCBI viral
        refseq database and RVDB.
    @param negativeTitleRegex: A regex that pathogen names must not match.
        Note that this matching is done on the final part of the protein title
        in square brackets, according to the convention used by the NCBI viral
        refseq database and RVDB.
    @param pathogenDataDir: The C{str} directory where per-pathogen information
        (e.g., collected reads across all samples) should be written. Will be
        created (in C{self.toHTML}) if it doesn't exist.
    @raise ValueError: If C{format_} is unknown.
    """

    VIRALZONE = "https://viralzone.expasy.org/search?query="
    ICTV = "https://talk.ictvonline.org/search-124283882/?q="
    READCOUNT_MARKER = "*READ-COUNT*"
    READ_AND_HSP_COUNT_STR_SEP = "/"

    def __init__(
        self,
        proteinGenomeDatabase,
        taxonomyDatabase,
        assetDir="out",
        sampleName=None,
        sampleNameRegex=None,
        format_="fasta",
        saveReadLengths=False,
        titleRegex=None,
        negativeTitleRegex=None,
        pathogenDataDir="pathogen-data",
    ):
        self._db = proteinGenomeDatabase
        self._taxdb = taxonomyDatabase
        self._assetDir = assetDir
        self._sampleName = sampleName
        self._sampleNameRegex = re.compile(sampleNameRegex) if sampleNameRegex else None
        if format_ in ("fasta", "fastq"):
            self._format = format_
        else:
            raise ValueError("format_ must be either 'fasta' or 'fastq'.")
        self._saveReadLengths = saveReadLengths

        if titleRegex or negativeTitleRegex:
            self.titleFilter = TitleFilter(
                positiveRegex=titleRegex, negativeRegex=negativeTitleRegex
            )
        else:
            self.titleFilter = None

        self._pathogenDataDir = pathogenDataDir

        # genomeAccessions will be a dict of dicts of dicts. The first
        # two keys will be a pathogen accession and a sample name. The
        # final dict will contain 'proteins' (a list of dicts) and
        # 'uniqueReadCount' (an int).
        self.genomeAccessions = defaultdict(dict)
        # sampleNames is keyed by sample name and will have values that hold
        # the sample's alignment panel index.html file.
        self.sampleNames = {}
        self.pathogenSampleFiles = PathogenSampleFiles(self, format_=format_)

    def _title(self, pathogenType):
        """
        Create a title summarizing the pathogens and samples.

        @param pathogenType: A C{str}, either 'viral', 'bacterial' or
            'generic'.
        @return: A C{str} title.
        """

        assert pathogenType in ("bacterial", "viral", "generic")

        nPathogens = len(self.genomeAccessions)
        nSamples = len(self.sampleNames)

        if pathogenType == "bacterial":
            what = "bacterium" if nPathogens == 1 else "bacteria"
        elif pathogenType == "viral":
            what = "virus%s" % ("" if nPathogens == 1 else "es")
        else:
            what = "pathogen%s" % ("" if nPathogens == 1 else "es")

        return "Proteins from %d %s were found in %d sample%s." % (
            nPathogens,
            what,
            nSamples,
            "" if nSamples == 1 else "s",
        )

    def addFile(self, filename, fp):
        """
        Read and record protein information for a sample.

        @param filename: A C{str} file name.
        @param fp: An open file pointer to read the file's data from.
        @raise ValueError: If information for a pathogen/protein/sample
            combination is given more than once.
        """
        if self._sampleName:
            sampleName = self._sampleName
        elif self._sampleNameRegex:
            match = self._sampleNameRegex.search(filename)
            if match:
                sampleName = match.group(1)
            else:
                sampleName = filename
        else:
            sampleName = filename

        outDir = join(dirname(filename), self._assetDir)

        self.sampleNames[sampleName] = join(outDir, "index.html")

        for index, proteinLine in enumerate(fp):
            (
                coverage,
                medianScore,
                bestScore,
                readCount,
                hspCount,
                proteinLength,
                longName,
            ) = proteinLine.split(None, 6)

            proteinInfo = self._db.findProtein(longName)
            if proteinInfo is None:
                try:
                    accession = self._db.proteinAccession(longName)
                except IndexError:
                    accession = longName
                # We could arguably just emit a warning here. This situation
                # arises (at least) when we are re-processing output from an
                # earlier run that used a different genome/protein
                # database. For example, the host specificity information about
                # a virus might change or the NCBI might withdraw it, causing
                # it to be excluded from a new database that we make. If an
                # now-not-present accession number appears in the DIAMOND or
                # alignment panel summary protein file, it will trigger this
                # error.
                #
                # For now I (Terry) have decided to keep things strict here and
                # raise an Exception. Otherwise I don't think there's any
                # guarantee that a warning to stderr would be seen, and only
                # issuing a warning would risk silently being in a situation
                # where nothing at all matched, e.g., due to passing an
                # incorrect database name. This error happens infrequently and
                # IMO it's better that we cause an error, force the user
                # (usually me, unfortunately) to investigate, clean up
                # properly, and re-run.
                raise ValueError(
                    "Could not find protein info for accession number %r "
                    "(extracted from %r). In the past, this hard-to-debug "
                    "(hence this long message!) error has resulted from using "
                    "a new genome/protein database to process results that "
                    "were generated based on an earlier version of the "
                    "database, in which case proteins that were present then "
                    "are not now in the database." % (accession, longName)
                )
            proteinName = proteinInfo["product"] or proteinInfo["gene"] or "unknown"
            proteinAccession = proteinInfo["accession"]

            genomeInfo = self._db.findGenome(longName)
            genomeName = genomeInfo["name"]
            genomeAccession = genomeInfo["accession"]

            # Ignore genomes with names we don't want.
            if (
                self.titleFilter
                and self.titleFilter.accept(genomeName) == TitleFilter.REJECT
            ):
                continue

            if sampleName not in self.genomeAccessions[genomeAccession]:
                self.genomeAccessions[genomeAccession][sampleName] = {
                    "proteins": {},
                    "uniqueReadCount": None,
                }

            proteins = self.genomeAccessions[genomeAccession][sampleName]["proteins"]

            # We should only receive one line of information for a given
            # genome/sample/protein combination.
            if proteinAccession in proteins:
                raise ValueError(
                    "Protein %r already seen for genome %r (%s) sample %r."
                    % (proteinAccession, genomeName, genomeAccession, sampleName)
                )

            readsFilename = join(outDir, "%s.%s" % (proteinAccession, self._format))

            if longName.startswith(
                SqliteIndexWriter.SEQUENCE_ID_PREFIX
                + SqliteIndexWriter.SEQUENCE_ID_SEPARATOR
            ):
                proteinURL = NCBISequenceLinkURL(longName, field=2)
                genomeURL = NCBISequenceLinkURL(longName, field=4)
            else:
                proteinURL = genomeURL = None

            proteinInfo = proteins[proteinAccession] = {
                "accession": proteinAccession,
                "bestScore": float(bestScore),
                "bluePlotFilename": join(outDir, "%s.png" % proteinAccession),
                "coverage": float(coverage),
                "readsFilename": readsFilename,
                "hspCount": int(hspCount),
                "index": index,
                "medianScore": float(medianScore),
                "outDir": outDir,
                "proteinLength": int(proteinLength),
                "proteinName": proteinName,
                "proteinURL": proteinURL,
                "genomeURL": genomeURL,
                "readCount": int(readCount),
            }

            if proteinInfo["readCount"] == proteinInfo["hspCount"]:
                proteinInfo["readAndHspCountStr"] = readCount
            else:
                proteinInfo["readAndHspCountStr"] = "%s%s%s" % (
                    readCount,
                    self.READ_AND_HSP_COUNT_STR_SEP,
                    hspCount,
                )

            if self._saveReadLengths:
                readsClass = FastaReads if self._format == "fasta" else FastqReads
                proteins[proteinName]["readLengths"] = tuple(
                    len(read) for read in readsClass(readsFilename)
                )

    def _computeUniqueReadCounts(self):
        """
        Add all pathogen / sample combinations to self.pathogenSampleFiles.

        This will make all de-duplicated (by id) FASTA/FASTQ files and store
        the number of de-duplicated reads into C{self.genomeAccessions}.
        """
        for genomeAccession, samples in self.genomeAccessions.items():
            for sampleName in samples:
                self.pathogenSampleFiles.add(genomeAccession, sampleName)

    def toStr(self, title=None, preamble=None, pathogenType="viral"):
        """
        Produce a string representation of the pathogen summary.

        @param title: The C{str} title for the output.
        @param preamble: The C{str} descriptive preamble, or C{None} if no
            preamble is needed.
        @param pathogenType: A C{str}, either 'viral', 'bacterial' or
            'generic'.

        @return: A C{str} suitable for printing.
        """
        # Note that the string representation contains much less
        # information than the HTML summary. E.g., it does not contain the
        # unique (de-duplicated, by id) read count, since that is only computed
        # when we are making combined FASTA files of reads matching a
        # pathogen.

        assert pathogenType in ("viral", "bacterial", "generic")

        title = title or "Summary of %s." % (
            "bacteria"
            if pathogenType == "bacterial"
            else ("viruses" if pathogenType == "viral" else "pathogens")
        )

        readCountGetter = itemgetter("readCount")
        result = []
        append = result.append

        result.extend((title, ""))
        if preamble:
            result.extend((preamble, ""))
        result.extend((self._title(pathogenType), ""))

        for genomeAccession, samples in self.genomeAccessions.items():
            genomeInfo = self._db.findGenome(genomeAccession)
            genomeName = genomeInfo["name"]
            sampleCount = len(samples)
            append(
                "%s (in %d sample%s)"
                % (genomeName, sampleCount, "" if sampleCount == 1 else "s")
            )
            for sampleName in sorted(samples):
                proteins = samples[sampleName]["proteins"]
                proteinCount = len(proteins)
                totalReads = sum(readCountGetter(p) for p in proteins.values())
                append(
                    "  %s (%d protein%s, %d read%s)"
                    % (
                        sampleName,
                        proteinCount,
                        "" if proteinCount == 1 else "s",
                        totalReads,
                        "" if totalReads == 1 else "s",
                    )
                )
                for proteinName in sorted(proteins):
                    append(
                        "    %(coverage).2f\t%(medianScore).2f\t"
                        "%(bestScore).2f\t%(readAndHspCountStr)3s\t"
                        "%(proteinName)s" % proteins[proteinName]
                    )
            append("")

        return "\n".join(result)

    def _genomeName(self, genomeAccession):
        """
        Get the name of a genome, given its accession number.

        @param genomeAccession: A C{str} pathogen accession number.
        @return: A C{str} genome name.
        """
        return self._db.findGenome(genomeAccession)["organism"]

    def _makeSampleSorter(self):
        """
        Make a function to sort sample names with, using the 3rd
        underscore-separated field of each name as an integer, if possible.
        """
        # Note: we could do this without the allSampleNamesHaveIntThirdField
        # variable by defining a function in the 'except' clause and adding an
        # 'else' to the 'for' loop, but that causes flake8 to complain that the
        # unused _key function (in the except) has been redefined (in the
        # else).
        allSampleNamesHaveIntThirdField = True
        for sampleName in self.sampleNames:
            try:
                int(sampleName.split("_", maxsplit=3)[2])
            except (IndexError, ValueError):
                allSampleNamesHaveIntThirdField = False
                break

        if allSampleNamesHaveIntThirdField:

            def _key(sampleName):
                return int(sampleName.split("_", maxsplit=3)[2])

        else:

            def _key(sampleName):
                return sampleName

        self.sampleSort = partial(sorted, key=_key)

    def toHTML(
        self,
        pathogenPanelFilename=None,
        readCountColors=None,
        minProteinFraction=0.0,
        minProteinCount=0,
        pathogenType="viral",
        title=None,
        preamble=None,
        sampleIndexFilename=None,
        omitVirusLinks=False,
        bootstrapTreeviewDir=None,
    ):
        """
        Produce an HTML string representation of the pathogen summary.

        @param pathogenPanelFilename: If not C{None}, a C{str} filename to
            write a pathogen panel PNG image to.
        @param readCountColors: Either a C{dark.colors.colorsForCounts}
            instance or C{None} for no read count coloring.
        @param minProteinFraction: The C{float} minimum fraction of proteins
            in a pathogen that must be matched by a sample in order for that
            pathogen to be displayed for that sample.
        @param minProteinCount: The C{int} minimum number of proteins
            in a pathogen that must be matched by a sample in order for that
            pathogen to be displayed for that sample.
        @param pathogenType: A C{str} giving the type of the pathogen involved,
            either 'bacterial', 'viral', or 'generic'.
        @param title: The C{str} title for the HTML page or C{None} to get a
            default generic title depending on whether a viral or bacterial
            database was matched against.
        @param preamble: The C{str} descriptive preamble for the HTML page, or
            C{None} if no preamble is needed.
        @param sampleIndexFilename: A C{str} filename to write a sample index
            file to. Lines in the file will have an integer index, a space, and
            then the sample name.
        @param omitVirusLinks: If C{True}, links to ICTV and ViralZone will be
            omitted in output.
        @param bootstrapTreeviewDir: A C{str} giving the directory where the
            bootstrap-treeview JS and CSS files may be found. Or C{None} if no
            bootstrap-treeview output should be generated.
        @return: An HTML C{str} suitable for printing.
        """
        if pathogenType == "bacterial":
            singular, plural = "bacterium", "bacteria"
        elif pathogenType == "viral":
            singular, plural = "virus", "viruses"
        elif pathogenType == "generic":
            singular, plural = "pathogen", "pathogens"
        else:
            raise ValueError(
                "Unrecognized pathogenType argument: %r. Value must be either "
                "'bacterial', 'viral', or 'generic'." % pathogenType
            )

        if not exists(self._pathogenDataDir):
            os.mkdir(self._pathogenDataDir)

        title = title or "Summary of " + plural

        self._makeSampleSorter()
        self._computeUniqueReadCounts()

        if sampleIndexFilename:
            with open(sampleIndexFilename, "w") as fp:
                self.pathogenSampleFiles.writeSampleIndex(fp)

        # Figure out if we have to delete some pathogens because the number
        # or fraction of its proteins that we have matches for is too low.
        if minProteinFraction > 0.0 or minProteinCount > 0:
            toDelete = defaultdict(list)
            for genomeAccession in self.genomeAccessions:
                genomeInfo = self._db.findGenome(genomeAccession)
                pathogenProteinCount = genomeInfo["proteinCount"]
                assert pathogenProteinCount > 0
                for s in self.genomeAccessions[genomeAccession]:
                    sampleProteinCount = len(
                        self.genomeAccessions[genomeAccession][s]["proteins"]
                    )
                    if sampleProteinCount < minProteinCount:
                        toDelete[genomeAccession].append(s)
                    else:
                        sampleProteinFraction = (
                            sampleProteinCount / pathogenProteinCount
                        )
                        if sampleProteinFraction < minProteinFraction:
                            toDelete[genomeAccession].append(s)

            for genomeAccession, samples in toDelete.items():
                for sample in samples:
                    del self.genomeAccessions[genomeAccession][sample]

        genomeAccessions = sorted(
            (
                genomeAccession
                for genomeAccession in self.genomeAccessions
                if len(self.genomeAccessions[genomeAccession]) > 0
            ),
            key=self._genomeName,
        )
        nPathogenNames = len(genomeAccessions)

        sampleNames = self.sampleSort(self.sampleNames)

        # Be very careful with commas in the following! Long lines that
        # should be continued unbroken must not end with a comma.
        result = [
            "<html>",
            "<head>",
            "<title>",
            title,
            "</title>",
            '<meta charset="UTF-8">',
            '<link rel="stylesheet"',
            'href="https://stackpath.bootstrapcdn.com/bootstrap/'
            '3.4.1/css/bootstrap.min.css"',
            'integrity="sha384-HSMxcRTRxnN+Bdg0JdbxYKrThecOKuH5z'
            'CYotlSAcp1+c8xmyTe9GYg1l9a69psu"',
            'crossorigin="anonymous">',
        ]

        if bootstrapTreeviewDir:
            result.append(
                '<link rel="stylesheet" href="%s/bootstrap-treeview.min.css">'
                % bootstrapTreeviewDir
            )

        result.extend(
            [
                "</head>",
                "<body>",
                "<script",
                'src="https://code.jquery.com/jquery-3.4.1.min.js"',
                'integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo="',
                'crossorigin="anonymous"></script>',
                "<script",
                'src="https://stackpath.bootstrapcdn.com/bootstrap/'
                '3.4.1/js/bootstrap.min.js"',
                'integrity="sha384-aJ21OjlMXNL5UyIl/XNwTMqvzeRMZH2w8c5cRVpzpU8Y5b'
                'ApTppSuUkhZXN0VxHd"',
                'crossorigin="anonymous"></script>',
            ]
        )

        if bootstrapTreeviewDir:
            result.append(
                '<script src="%s/bootstrap-treeview.min.js"></script>'
                % bootstrapTreeviewDir
            )

        result.extend(
            [
                "<style>",
                """\
            body {
                margin-left: 2%;
                margin-right: 2%;
            }
            hr {
                display: block;
                margin-top: 0.5em;
                margin-bottom: 0.5em;
                margin-left: auto;
                margin-right: auto;
                border-style: inset;
                border-width: 1px;
            }
            p.pathogen {
                margin-top: 10px;
                margin-bottom: 3px;
            }
            p.sample {
                margin-top: 10px;
                margin-bottom: 3px;
            }
            .sample {
                margin-top: 5px;
                margin-bottom: 2px;
            }
            ul {
                margin-bottom: 2px;
            }
            .indented {
                margin-left: 2em;
            }
            .sample-name {
                font-size: 115%;
                font-weight: bold;
            }
            .pathogen-name {
                font-size: 115%;
                font-weight: bold;
            }
            .index-name {
                font-weight: bold;
            }
            .index {
                font-size: small;
            }
            .index-letter {
                font-size: 115%;
                font-weight: bold;
            }
            .host {
                font-size: small;
            }
            .taxonomy {
                font-size: small;
            }
            .protein-name {
            }
            .stats {
                font-family: "Courier New", Courier, monospace;
                white-space: pre;
            }
            .protein-list {
                margin-top: 2px;
            }""",
                "</style>",
                "</head>",
                "<body>",
            ]
        )

        append = result.append

        proteinFieldsDescription = self._help(readCountColors, result)

        append("<h2>%s</h2>" % title)
        append(self._title(pathogenType))
        if preamble:
            append(preamble)

        if minProteinFraction > 0.0:
            append("<p>")
            percent = minProteinFraction * 100.0
            if nPathogenNames < len(self.genomeAccessions):
                if nPathogenNames == 1:
                    append(
                        "Pathogen protein fraction filtering has been "
                        "applied, so information on only 1 pathogen is "
                        "displayed. This is the only pathogen for which at "
                        "least one sample matches at least %.2f%% of the "
                        "pathogen proteins." % percent
                    )
                else:
                    append(
                        "Pathogen protein fraction filtering has been "
                        "applied, so information on only %d pathogens is "
                        "displayed. These are the only pathogens for which "
                        "at least one sample matches at least %.2f%% of "
                        "the pathogen proteins." % (nPathogenNames, percent)
                    )
            else:
                append(
                    "Pathogen protein fraction filtering was applied, "
                    "but all pathogens have at least %.2f%% of their "
                    "proteins matched by at least one sample." % percent
                )
            append("</p>")

        if pathogenPanelFilename and genomeAccessions:
            self.pathogenPanel(pathogenPanelFilename)
            append("<p>")
            append(
                '<a href="%s">Panel showing read count per pathogen, '
                "per sample.</a>" % pathogenPanelFilename
            )
            append(
                "Red vertical bars indicate samples with an unusually "
                "high read count."
            )
            append("</p>")

        result.extend(proteinFieldsDescription)

        append('<p style="margin-top: 10px;">Global: ')
        append(
            '<button type="button" class="btn btn-default btn-sm" '
            'id="expand-all-button">Expand all</button>'
        )

        append(
            '<button type="button" class="btn btn-default btn-sm" '
            'id="collapse-all-button">Collapse all</button>'
        )
        append("</p>")

        append(
            """
        <script>
        $("#expand-all-button").click(function(){
            $(".collapse").collapse("show");
        });
        $("#collapse-all-button").click(function(){
            $(".collapse").collapse("hide");
        });
        </script>
        """
        )

        append("<h2>Indices</h2>")
        self._sampleIndex(sampleNames, result)
        self._pathogenIndex(genomeAccessions, result, singular, plural)

        self._samplesToHTML(
            result,
            pathogenType,
            omitVirusLinks,
            sampleNames,
            readCountColors,
            singular,
            plural,
        )

        self._pathogensToHTML(
            result,
            pathogenType,
            genomeAccessions,
            omitVirusLinks,
            readCountColors,
            bootstrapTreeviewDir,
            plural,
        )

        append("</body>")
        append("</html>")

        return "\n".join(result)

    def _samplesToHTML(
        self,
        result,
        pathogenType,
        omitVirusLinks,
        sampleNames,
        readCountColors,
        singular,
        plural,
    ):
        """
        Write all samples (with pathogens (with proteins)).
        """
        append = result.append
        append("<h2>Samples</h2>")

        for sampleName in sampleNames:
            samplePathogenAccessions = sorted(
                (
                    accession
                    for accession in self.genomeAccessions
                    if sampleName in self.genomeAccessions[accession]
                ),
                key=self._genomeName,
            )

            append("<div>")

            append(
                '<button type="button" class="btn btn-default btn-sm" '
                'data-toggle="collapse" data-target="#sample-%s-collapse">'
                '<span class="glyphicon glyphicon-plus"></span></button>' % sampleName
            )

            if len(samplePathogenAccessions):
                append(
                    '<a id="sample-%s"></a>'
                    '<span class="sample"><span class="sample-name">%s</span> '
                    "matched proteins from %d %s, "
                    '<a href="%s">panel</a>.</span>'
                    % (
                        sampleName,
                        sampleName,
                        len(samplePathogenAccessions),
                        (singular if len(samplePathogenAccessions) == 1 else plural),
                        self.sampleNames[sampleName],
                    )
                )
            else:
                append(
                    '<a id="sample-%s"></a>'
                    '<span class="sample">'
                    '<span class="sample-name">%s</span> '
                    "did not match anything.</span>" % (sampleName, sampleName)
                )
                continue

            append("</div>")
            append('<div class="collapse" id="sample-%s-collapse">' % sampleName)

            for genomeAccession in samplePathogenAccessions:
                genomeInfo = self._db.findGenome(genomeAccession)
                readsFileName = self.pathogenSampleFiles.lookup(
                    genomeAccession, sampleName
                )
                proteins = self.genomeAccessions[genomeAccession][sampleName][
                    "proteins"
                ]
                uniqueReadCount = self.genomeAccessions[genomeAccession][sampleName][
                    "uniqueReadCount"
                ]
                proteinCount = len(proteins)
                pathogenProteinCount = genomeInfo["proteinCount"]
                proteinCountStr = "%d/%d protein%s" % (
                    proteinCount,
                    pathogenProteinCount,
                    "" if pathogenProteinCount == 1 else "s",
                )

                pathogenLinksHTML = " (%s" % NCBISequenceLink(genomeAccession)

                if pathogenType == "viral" and not omitVirusLinks:
                    quoted = quote(genomeInfo["organism"])
                    pathogenLinksHTML += (
                        ', <a href="%s%s">ICTV</a>, ' '<a href="%s%s">ViralZone</a>)'
                    ) % (self.ICTV, quoted, self.VIRALZONE, quoted)
                else:
                    pathogenLinksHTML += ")"

                append(
                    '<p class="sample indented">'
                    '<a href="#pathogen-%s">%s</a> %s %s, '
                    '<a href="%s">%d read%s</a>:</p>'
                    % (
                        genomeAccession,
                        genomeInfo["organism"],
                        pathogenLinksHTML,
                        proteinCountStr,
                        readsFileName,
                        uniqueReadCount,
                        "" if uniqueReadCount == 1 else "s",
                    )
                )
                append('<ul class="protein-list indented">')
                for proteinAccession in sorted(proteins):
                    proteinMatch = proteins[proteinAccession]
                    append(
                        "<li>"
                        '<span class="stats">'
                        "%(coverage).2f %(medianScore)6.2f %(bestScore)6.2f "
                        % proteinMatch
                    )

                    self._appendNoSpace(
                        readCountText(
                            readCountColors,
                            proteinMatch["readCount"],
                            f"{proteinMatch['readAndHspCountStr']:4s}",
                        ),
                        result,
                    )

                    self._appendNoSpace(
                        "</span> "
                        '<span class="protein-name">'
                        "%(proteinName)s"
                        "</span> "
                        "(%(proteinLength)d aa," % proteinMatch,
                        result,
                    )

                    if proteinMatch["proteinURL"]:
                        append(
                            '<a href="%s">%s</a>, '
                            % (proteinMatch["proteinURL"], proteinMatch["accession"])
                        )

                    append(
                        '<a href="%(bluePlotFilename)s">blue plot</a>, '
                        '<a href="%(readsFilename)s">reads</a>)' % proteinMatch
                    )

                    append("</li>")

                append("</ul>")
            append("</div>")

    def _pathogensToHTML(
        self,
        result,
        pathogenType,
        genomeAccessions,
        omitVirusLinks,
        readCountColors,
        bootstrapTreeviewDir,
        plural,
    ):
        """
        Write all pathogens (with samples (with proteins)).
        """
        append = result.append
        append("<h2>%s</h2>" % plural.title())

        if bootstrapTreeviewDir:
            # A <div> to hold the taxonomy tree.
            append('<div id="tree"></div>')

        taxonomyHierarchy = Hierarchy()

        for genomeAccession in genomeAccessions:
            samples = self.genomeAccessions[genomeAccession]
            sampleCount = len(samples)
            genomeInfo = self._db.findGenome(genomeAccession)
            pathogenProteinCount = genomeInfo["proteinCount"]

            lineage = (
                None
                if genomeInfo["taxonomyId"] is None
                else self._taxdb.lineage(genomeInfo["taxonomyId"])
            )

            if lineage:
                taxonomyHierarchy.add(lineage, genomeAccession)
                lineageHTML = ", ".join(lineageTaxonomyLinks(lineage))
            else:
                lineageHTML = ""

            pathogenLinksHTML = " %s, %s" % (
                genomeInfo["databaseName"],
                NCBISequenceLink(genomeAccession),
            )

            if pathogenType == "viral" and not omitVirusLinks:
                quoted = quote(genomeInfo["organism"])
                pathogenLinksHTML += (
                    ', <a href="%s%s">ICTV</a>, <a href="%s%s">ViralZone</a>.'
                ) % (self.ICTV, quoted, self.VIRALZONE, quoted)
            else:
                pathogenLinksHTML += "."

            proteinCountStr = " %d protein%s" % (
                pathogenProteinCount,
                "" if pathogenProteinCount == 1 else "s",
            )

            pathogenReadsFilename = join(
                self._pathogenDataDir,
                "pathogen-%s.%s" % (genomeAccession, self._format),
            )

            pathogenReadsFp = open(pathogenReadsFilename, "w")
            pathogenReadCount = 0

            append("<div>")  # Button and following summary.

            append(
                '<button type="button" class="btn btn-default btn-sm" '
                'data-toggle="collapse" '
                'data-target="#pathogen-%s-collapse">'
                '<span class="glyphicon glyphicon-plus"></span></button>'
                % genomeAccession.replace(".", "-")
            )

            append(
                '<a id="pathogen-%s"></a>'
                '<span class="pathogen">'
                '<span class="pathogen-name">%s</span> '
                '<span class="host">(%s)</span>'
                "<br/>%d nt, %s, "
                "matched by %d sample%s, "
                '%s <a href="%s">reads</a> in total. '
                "%s"
                '<br/><span class="taxonomy">Taxonomy: %s.</span>'
                "</span>"
                % (
                    genomeAccession,
                    genomeInfo["organism"],
                    genomeInfo.get("host") or "unknown host",
                    genomeInfo["length"],
                    proteinCountStr,
                    sampleCount,
                    "" if sampleCount == 1 else "s",
                    self.READCOUNT_MARKER,
                    pathogenReadsFilename,
                    pathogenLinksHTML,
                    lineageHTML,
                )
            )

            # Remember where we are in the output result so we can fill in
            # the total read count once we have processed all samples for
            # this pathogen. Not nice, I know.
            pathogenReadCountLineIndex = len(result) - 1

            append("</div>")  # End of button summary.

            append(
                '<div class="collapse" id="pathogen-%s-collapse">'
                % genomeAccession.replace(".", "-")
            )

            for sampleName in self.sampleSort(samples):
                readsFileName = self.pathogenSampleFiles.lookup(
                    genomeAccession, sampleName
                )

                # Copy the read data from the per-sample reads for this
                # pathogen into the per-pathogen file of reads.
                with open(readsFileName) as readsFp:
                    while True:
                        data = readsFp.read(4096)
                        if data:
                            pathogenReadsFp.write(data)
                        else:
                            break

                proteins = samples[sampleName]["proteins"]
                proteinCount = len(proteins)
                uniqueReadCount = samples[sampleName]["uniqueReadCount"]
                pathogenReadCount += uniqueReadCount
                proteinCountHTML = "%d protein%s, " % (
                    proteinCount,
                    "" if proteinCount == 1 else "s",
                )

                append(
                    '<p class="sample indented">'
                    'Sample <a href="#sample-%s">%s</a> '
                    '(%s<a href="%s">%d '
                    'read%s</a>, <a href="%s">panel</a>).</p>'
                    % (
                        sampleName,
                        sampleName,
                        proteinCountHTML,
                        readsFileName,
                        uniqueReadCount,
                        "" if uniqueReadCount == 1 else "s",
                        self.sampleNames[sampleName],
                    )
                )

                append('<ul class="protein-list indented">')
                for proteinName in sorted(proteins):
                    proteinMatch = proteins[proteinName]
                    append(
                        "<li>"
                        '<span class="stats">'
                        "%(coverage).2f %(medianScore)6.2f %(bestScore)6.2f "
                        % proteinMatch
                    )

                    self._appendNoSpace(
                        readCountText(
                            readCountColors,
                            proteinMatch["readCount"],
                            f"{proteinMatch['readAndHspCountStr']:4s}",
                        ),
                        result,
                    )

                    if self._saveReadLengths:
                        self._appendNoSpace(
                            " (%s)"
                            % ", ".join(map(str, sorted(proteinMatch["readLengths"]))),
                            result,
                        )

                    self._appendNoSpace(
                        "</span> "
                        '<span class="protein-name">'
                        "%(proteinName)s"
                        "</span> "
                        "(%(proteinLength)d aa," % proteinMatch,
                        result,
                    )

                    if proteinMatch["proteinURL"]:
                        append(
                            '<a href="%s">%s</a>, '
                            % (proteinMatch["proteinURL"], proteinMatch["accession"])
                        )

                    append(
                        '<a href="%(bluePlotFilename)s">blue plot</a>, '
                        '<a href="%(readsFilename)s">reads</a>)' % proteinMatch
                    )

                    append("</li>")

                append("</ul>")

            append("</div>")

            pathogenReadsFp.close()

            # Sanity check there's a read count marker text in our output
            # where we expect it.
            readCountLine = result[pathogenReadCountLineIndex]
            if readCountLine.find(self.READCOUNT_MARKER) == -1:
                raise ValueError(
                    "Could not find pathogen read count marker (%s) in result "
                    "index %d text (%s)."
                    % (self.READCOUNT_MARKER, pathogenReadCountLineIndex, readCountLine)
                )

            # Put the read count into the pathogen summary line we wrote
            # earlier, replacing the read count marker with the correct text.
            readCountLine = readCountLine.replace(
                self.READCOUNT_MARKER, readCountText(readCountColors, pathogenReadCount)
            )
            if pathogenReadCount == 1:
                # Horrible hack to make 'reads' be singular if there's only 1.
                readCountLine = readCountLine.replace(
                    '">reads</a> in total.', '">read</a> in total.'
                )

            result[pathogenReadCountLineIndex] = readCountLine

        if bootstrapTreeviewDir:
            append(
                """
                <script>
                $(document).ready(function(){
                    var tree = %s;
                    $('#tree').treeview({
                        data: tree,
                        enableLinks: true,
                        levels: 0,
                    });
                });
                </script>
            """
                % taxonomyHierarchy.toJSON()
            )

    def _help(self, readCountColors, result):
        append = result.append

        append(
            """
        <script>
            $(document).ready(function(){
                $("#help-button").click(function(){
                    var self=$(this);
                    if (self.val() === "Show"){
                        self.val("Hide");
                    }
                    else {
                        self.val("Show");
                    }
                    $("#help-details").toggle();
                });
            });
        </script>
        """
        )

        if readCountColors:
            levels = []
            append("<style>")
            for threshold, color in readCountColors.colors:
                klass = readCountColors.thresholdToCssName(threshold)
                append(".%s { color: %s; font-weight: bold; }" % (klass, color))
                levels.append('<span class="%s">%d</span>' % (klass, threshold))
            append("</style>")
            readCountColorLegend = " Color levels: " + ", ".join(reversed(levels)) + "."
        else:
            readCountColorLegend = ""

        proteinFieldsDescription = [
            'Help: <button type="button" class="btn btn-default btn-sm" '
            'id="help-button">Show</button><br>',
            '<div id="help-details" style="display:none;">',
            "In all bullet point protein lists below, there are the following "
            "numeric fields:",
            "<ol>",
            "<li>Coverage fraction.</li>",
            "<li>Median bit score.</li>",
            "<li>Best bit score.</li>",
            "<li>Read count (if read and HSP counts differ, ",
            (
                'both are given, separated by "%s").%s</li>'
                % (self.READ_AND_HSP_COUNT_STR_SEP, readCountColorLegend)
            ),
        ]

        if self._saveReadLengths:
            proteinFieldsDescription.append(
                "<li>All read lengths (in parentheses).</li>"
            )

        proteinFieldsDescription.extend(
            [
                "</ol>",
                "</div>",
            ]
        )

        return proteinFieldsDescription

    def _appendNoSpace(self, s, result):
        assert result, "Cannot append %r to empty result list" % s
        result[-1] += s

    def _sampleIndex(self, sampleNames, result):
        """
        Write a linked table of contents by sample.
        """
        append = result.append

        if len(sampleNames) == 1:
            title = "Sample"
        else:
            title = "Samples (%d)" % len(sampleNames)

        append('<p><span class="index-name">%s:</span>' % title)
        append('<span class="index">')
        for count, sampleName in enumerate(sampleNames, start=1):
            append(
                '<span class="index-letter">%d</span> '
                '<a href="#sample-%s">%s</a>' % (count, sampleName, sampleName)
            )
            append("&middot;")
        # Get rid of final middle dot and add a period.
        result.pop()
        self._appendNoSpace(".", result)
        append("</span></p>")

    def _pathogenIndex(self, genomeAccessions, result, singular, plural):
        """
        Create a linked table of contents by pathogen.
        """
        append = result.append

        if len(genomeAccessions) == 1:
            title = singular.title()
        else:
            title = "%s (%d)" % (plural.title(), len(genomeAccessions))

        append('<p><span class="index-name">%s:</span>' % title)
        append('<span class="index">')
        lastLetter = None
        for genomeAccession in genomeAccessions:
            genomeInfo = self._db.findGenome(genomeAccession)
            organism = genomeInfo["organism"]
            letter = organism[0]
            if letter != lastLetter:
                append('<span class="index-letter">%s</span>' % letter)
                lastLetter = letter
            append(
                '<a href="#pathogen-%s">%s</a>'
                % (genomeAccession, genomeInfo["organism"])
            )
            append("&middot;")
        # Get rid of final middle dot and add a period.
        result.pop()
        self._appendNoSpace(".", result)
        append("</span></p>")

    def _pathogenSamplePlot(self, genomeAccession, sampleNames, ax):
        """
        Make an image of a graph giving pathogen read count (Y axis) versus
        sample id (X axis).

        @param genomeAccession: A C{str} pathogen accession number.
        @param sampleNames: A sorted C{list} of sample names.
        @param ax: A matplotlib C{axes} instance.
        """
        readCounts = []
        for sampleName in sampleNames:
            try:
                readCount = self.genomeAccessions[genomeAccession][sampleName][
                    "uniqueReadCount"
                ]
            except KeyError:
                readCount = 0
            readCounts.append(readCount)

        highlight = "r"
        normal = "gray"
        sdMultiple = 2.5
        minReadsForHighlighting = 10
        highlighted = []

        if len(readCounts) == 1:
            if readCounts[0] > minReadsForHighlighting:
                color = [highlight]
                highlighted.append(sampleNames[0])
            else:
                color = [normal]
        else:
            mean = np.mean(readCounts)
            sd = np.std(readCounts)
            color = []
            for readCount, sampleName in zip(readCounts, sampleNames):
                if (
                    readCount > (sdMultiple * sd) + mean
                    and readCount >= minReadsForHighlighting
                ):
                    color.append(highlight)
                    highlighted.append(sampleName)
                else:
                    color.append(normal)

        nSamples = len(sampleNames)
        x = np.arange(nSamples)
        yMin = np.zeros(nSamples)
        ax.set_xticks([])
        ax.set_xlim((-0.5, nSamples - 0.5))
        ax.vlines(x, yMin, readCounts, color=color)

        genomeInfo = self._db.findGenome(genomeAccession)

        if highlighted:
            title = "%s\nIn red: %s" % (
                genomeInfo["organism"],
                fill(", ".join(highlighted), 50),
            )
        else:
            # Add a newline to keep the first line of each title at the
            # same place as those titles that have an "In red:" second
            # line.
            title = genomeInfo["organism"] + "\n"

        ax.set_title(title, fontsize=10)
        ax.tick_params(axis="both", which="major", labelsize=8)
        ax.tick_params(axis="both", which="minor", labelsize=6)

    def pathogenPanel(self, filename):
        """
        Make a panel of images, with each image being a graph giving pathogen
        de-duplicated (by id) read count (Y axis) versus sample id (X axis).

        @param filename: A C{str} file name to write the image to.
        """
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        self._computeUniqueReadCounts()
        genomeAccessions = sorted(self.genomeAccessions)
        sampleNames = sorted(self.sampleNames)

        cols = 5
        rows = int(len(genomeAccessions) / cols) + (
            0 if len(genomeAccessions) % cols == 0 else 1
        )
        figure, ax = plt.subplots(rows, cols, squeeze=False)

        coords = dimensionalIterator((rows, cols))

        for genomeAccession in genomeAccessions:
            row, col = next(coords)
            self._pathogenSamplePlot(genomeAccession, sampleNames, ax[row][col])

        # Hide the final panel graphs (if any) that have no content. We do
        # this because the panel is a rectangular grid and some of the
        # plots at the end of the last row may be unused.
        for row, col in coords:
            ax[row][col].axis("off")

        figure.suptitle(
            "Per-sample read count for %d pathogen%s and %d sample%s.\n\n"
            % (
                len(genomeAccessions),
                "" if len(genomeAccessions) == 1 else "s",
                len(sampleNames),
                "" if len(sampleNames) == 1 else "s",
            ),
            fontsize=18,
        )
        figure.set_size_inches(5.0 * cols, 2.0 * rows, forward=True)
        plt.subplots_adjust(hspace=0.4)

        try:
            figure.savefig(filename)
        except ValueError as e:
            print(
                f"WARNING! Could not save pathogens panel figure: {str(e)!r}. "
                "That file has not been created and therefore the link to it "
                "from the results HTML will be broken.",
                file=sys.stderr,
            )


class _Genome:
    """
    Hold genome information, mirroring the attributes of a BioPython
    GenBank record.

    @param d: A C{dict} holding genome information (see below).
    """

    def __init__(self, d):
        self.id = d["id"]
        self.description = d["name"]
        self.seq = d["sequence"]
        self.annotations = {}
        self.lineage = [LineageElement(*lineage) for lineage in d.get("lineage", [])]
        self.features = [_GenomeFeature(f) for f in d["features"]]


class _GenomeLocation:
    """
    Hold genome feature location information, mirroring the attributes of a
    BioPython GenBank record.

    @param start: An C{int} start location.
    @param end: An C{int} stop location.
    @param strand: The C{int} strand, either 1 for forward or 0 for reverse.
    """

    def __init__(self, start, end, strand):
        self.start = start
        self.end = end
        self.strand = strand

    def __str__(self):
        return "[%d:%d](%s)" % (self.start, self.end, "+" if self.strand == 1 else "-")


class _GenomeFeature:
    """
    Hold genome feature information, mirroring the attributes of a BioPython
    GenBank record.

    @param d: A C{dict} holding genome feature information.
    """

    def __init__(self, d):
        self.type = d["type"]
        self.qualifiers = d["qualifiers"]
        self.strand = 1
        location = d["qualifiers"]["location"]
        self.location = _GenomeLocation(
            location["start"], location["stop"], self.strand
        )


class SqliteIndexWriter:
    """
    Create or update an Sqlite3 database holding information about proteins and
    the genomes they come from.

    @param dbFilename: A C{str} file name containing an sqlite3 database. If
        the file does not exist it will be created. The special string
        ":memory:" can be used to create an in-memory database.
    @param fastaFp: A file-pointer to which the protein FASTA is written.
    """

    PROTEIN_ACCESSION_FIELD = 2
    GENOME_ACCESSION_FIELD = 4
    SEQUENCE_ID_PREFIX = "civ"
    SEQUENCE_ID_SEPARATOR = "|"

    def __init__(self, dbFilename, fastaFp=sys.stdout):
        self._connection = sqliteConnect(dbFilename)
        self._fastaFp = fastaFp

        cur = self._connection.cursor()
        cur.executescript(
            """
            CREATE TABLE IF NOT EXISTS proteins (
                accession VARCHAR UNIQUE PRIMARY KEY,
                genomeAccession VARCHAR NOT NULL,
                sequence VARCHAR NOT NULL,
                length INTEGER NOT NULL,
                offsets VARCHAR NOT NULL,
                forward INTEGER NOT NULL,
                circular INTEGER NOT NULL,
                rangeCount INTEGER NOT NULL,
                gene VARCHAR,
                note VARCHAR,
                product VARCHAR,
                FOREIGN KEY (genomeAccession)
                    REFERENCES genomes (accession)
            );

            CREATE TABLE IF NOT EXISTS genomes (
                accession VARCHAR UNIQUE PRIMARY KEY,
                organism VARCHAR NOT NULL,
                name VARCHAR NOT NULL,
                sequence VARCHAR NOT NULL,
                length INTEGER NOT NULL,
                proteinCount INTEGER NOT NULL,
                host VARCHAR,
                note VARCHAR,
                taxonomyId INTEGER,
                databaseName VARCHAR
            );
            """
        )
        self._connection.commit()

    def addGenBankFile(
        self,
        filename,
        taxonomyDatabase,
        dnaOnly=False,
        rnaOnly=False,
        allowedTaxonomicRanks=None,
        minGenomeLength=None,
        maxGenomeLength=None,
        excludeExclusiveHosts=None,
        excludeFungusOnlyViruses=False,
        excludePlantOnlyViruses=False,
        databaseName=None,
        proteinSource="GENBANK",
        genomeSource="GENBANK",
        duplicationPolicy="error",
        logfp=None,
    ):
        """
        Add proteins from a GenBank file.

        @param filename: A C{str} file name, with the file in GenBank format
            (see https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html).
        @param taxonomyDatabase: A taxonomy database. Must be given if
            C{dnaOnly} is C{True} or C{rnaOnly} is C{True} or
            C{excludeExclusiveHosts} is not C{None}.
        @param dnaOnly: If C{True}, only include DNA viruses.
        @param rnaOnly: If C{True}, only include RNA viruses.
        @param allowedTaxonomicRanks: If not C{None}, a set of case-insensitive
            name, rank C{str} 2-tuples. E.g.,
                set((
                    ("Nidovirales", "order"),
                    ("Retroviridae", "family"),
                ))
            Viruses will be included only if they match at least one of
            the name, rank pairs.
        @param minGenomeLength: If not C{None}, genomes of a length shorter
            than this should not be added.
        @param maxGenomeLength: If not C{None}, genomes of a length greater
            than this should not be added.
        @param excludeExclusiveHosts: Either C{None} or a set of host types
            that should cause a genome to be excluded if the genome only
            has a single host and it is in C{excludeExclusiveHosts}.
        @param excludeFungusOnlyViruses: If C{True}, do not include fungus-only
            viruses.
        @param excludePlantOnlyViruses: If C{True}, do not include plant-only
            viruses.
        @param databaseName: A C{str} indicating the database the records
            in C{filename} came from (e.g., 'refseq' or 'RVDB').
        @param proteinSource: A C{str} giving the source of the protein
            accession number. This becomes part of the sequence id printed
            in the protein FASTA output.
        @param genomeSource: A C{str} giving the source of the genome
            accession number. This becomes part of the sequence id printed
            in the protein FASTA output.
        @param duplicationPolicy: A C{str} indicating what to do if a
            to-be-inserted accession number is already present in the database.
            "error" results in a ValueError being raised, "ignore" means ignore
            the duplicate. It should also be possible to update (i.e., replace)
            but that is not supported yet.
        @param logfp: If not C{None}, a file pointer to write verbose
            progress output to.
        @raise DatabaseDuplicationError: If a duplicate accession number is
            encountered and C{duplicationPolicy} is 'error'.
        @return: A tuple containing two C{int}s: the number of genome sequences
            in the added file and the total number of proteins found.
        """

        def lineageFetcher(genome):
            return taxonomyDatabase.lineage(genome.id)

        with asHandle(filename) as fp:
            with self._connection:
                genomes = SeqIO.parse(fp, "gb")
                return self._addGenomes(
                    genomes,
                    taxonomyDatabase,
                    lineageFetcher,
                    dnaOnly=dnaOnly,
                    rnaOnly=rnaOnly,
                    allowedTaxonomicRanks=allowedTaxonomicRanks,
                    minGenomeLength=minGenomeLength,
                    maxGenomeLength=maxGenomeLength,
                    excludeExclusiveHosts=excludeExclusiveHosts,
                    excludeFungusOnlyViruses=excludeFungusOnlyViruses,
                    excludePlantOnlyViruses=excludePlantOnlyViruses,
                    databaseName=databaseName,
                    proteinSource=proteinSource,
                    genomeSource=genomeSource,
                    duplicationPolicy=duplicationPolicy,
                    logfp=logfp,
                )

    def addJSONFile(
        self,
        filename,
        taxonomyDatabase,
        dnaOnly=False,
        rnaOnly=False,
        allowedTaxonomicRanks=None,
        minGenomeLength=None,
        maxGenomeLength=None,
        excludeExclusiveHosts=None,
        excludeFungusOnlyViruses=False,
        excludePlantOnlyViruses=False,
        databaseName=None,
        proteinSource="GENBANK",
        genomeSource="GENBANK",
        duplicationPolicy="error",
        logfp=None,
    ):
        """
        Add proteins from a JSON infor file.

        @param filename: A C{str} file name, in JSON format.
        @param taxonomyDatabase: A taxonomy database. Must be given if
            C{dnaOnly} is C{True} or C{rnaOnly} is C{True} or
            C{excludeExclusiveHosts} is not C{None}.
        @param dnaOnly: If C{True}, only include DNA viruses.
        @param rnaOnly: If C{True}, only include RNA viruses.
        @param allowedTaxonomicRanks: If not C{None}, a set of case-insensitive
            name, rank C{str} 2-tuples. E.g.,
                set((
                    ("Nidovirales", "order"),
                    ("Retroviridae", "family"),
                ))
            Viruses will be included only if they match at least one of
            the name, rank pairs.
        @param minGenomeLength: If not C{None}, genomes of a length shorter
            than this should not be added.
        @param maxGenomeLength: If not C{None}, genomes of a length greater
            than this should not be added.
        @param excludeExclusiveHosts: Either C{None} or a set of host types
            that should cause a genome to be excluded if the genome only
            has a single host and it is in C{excludeExclusiveHosts}.
        @param excludeFungusOnlyViruses: If C{True}, do not include fungus-only
            viruses.
        @param excludePlantOnlyViruses: If C{True}, do not include plant-only
            viruses.
        @param databaseName: A C{str} indicating the database the records
            in C{filename} came from (e.g., 'refseq' or 'RVDB').
        @param proteinSource: A C{str} giving the source of the protein
            accession number. This becomes part of the sequence id printed
            in the protein FASTA output.
        @param genomeSource: A C{str} giving the source of the genome
            accession number. This becomes part of the sequence id printed
            in the protein FASTA output.
        @param duplicationPolicy: A C{str} indicating what to do if a
            to-be-inserted accession number is already present in the database.
            "error" results in a ValueError being raised, "ignore" means ignore
            the duplicate. It should also be possible to update (i.e., replace)
            but that is not supported yet.
        @param logfp: If not C{None}, a file pointer to write verbose
            progress output to.
        @raise DatabaseDuplicationError: If a duplicate accession number is
            encountered and C{duplicationPolicy} is 'error'.
        @return: A tuple containing two C{int}s: the number of genome sequences
            in the added file and the total number of proteins found.
        """

        def lineageFetcher(genome):
            return genome.lineage

        with asHandle(filename) as fp:
            genome = _Genome(load(fp))

        with self._connection:
            return self._addGenomes(
                [genome],
                taxonomyDatabase,
                lineageFetcher,
                dnaOnly=dnaOnly,
                rnaOnly=rnaOnly,
                allowedTaxonomicRanks=None,
                minGenomeLength=minGenomeLength,
                maxGenomeLength=maxGenomeLength,
                excludeExclusiveHosts=excludeExclusiveHosts,
                excludeFungusOnlyViruses=excludeFungusOnlyViruses,
                excludePlantOnlyViruses=excludePlantOnlyViruses,
                databaseName=databaseName,
                proteinSource=proteinSource,
                genomeSource=genomeSource,
                duplicationPolicy=duplicationPolicy,
                logfp=logfp,
            )

    def _addGenomes(
        self,
        genomes,
        taxonomyDatabase,
        lineageFetcher,
        dnaOnly=False,
        rnaOnly=False,
        allowedTaxonomicRanks=None,
        minGenomeLength=None,
        maxGenomeLength=None,
        excludeExclusiveHosts=None,
        excludeFungusOnlyViruses=False,
        excludePlantOnlyViruses=False,
        databaseName=None,
        proteinSource="GENBANK",
        genomeSource="GENBANK",
        duplicationPolicy="error",
        logfp=None,
    ):
        """
        Add a bunch of genomes.

        @param genomes: An iterable of genomes. These are either genomes
            returned by BioPython's GenBank parser or instances of C{_Genome}.
        @param taxonomyDatabase: A taxonomy database.
        @param lineageFetcher: A function that takes a genome and returns a
            C{tuple} of the taxonomic categories of the genome. Each
            tuple element is a 3-tuple of (C{int}, C{str}, C{str}) giving a
            taxonomy id a (scientific) name, and the rank (species, genus,
            etc). I.e., as returned by L{dark.taxonomy.LineageFetcher.lineage}.
        @param dnaOnly: If C{True}, only include DNA viruses.
        @param rnaOnly: If C{True}, only include RNA viruses.
        @param allowedTaxonomicRanks: If not C{None}, a set of case-insensitive
            name, rank C{str} 2-tuples. E.g.,
                set((
                    ("Nidovirales", "order"),
                    ("Retroviridae", "family"),
                ))
            Viruses will be included only if they match at least one of
            the name, rank pairs.
        @param minGenomeLength: If not C{None}, genomes of a length shorter
            than this should not be added.
        @param maxGenomeLength: If not C{None}, genomes of a length greater
            than this should not be added.
        @param excludeExclusiveHosts: Either C{None} or a set of host types
            that should cause a genome to be excluded if the genome only
            has a single host and it is in C{excludeExclusiveHosts}.
        @param excludeFungusOnlyViruses: If C{True}, do not include fungus-only
            viruses.
        @param excludePlantOnlyViruses: If C{True}, do not include plant-only
            viruses.
        @param databaseName: A C{str} indicating the database the records
            in C{filename} came from (e.g., 'refseq' or 'RVDB').
        @param proteinSource: A C{str} giving the source of the protein
            accession number. This becomes part of the sequence id printed
            in the protein FASTA output.
        @param genomeSource: A C{str} giving the source of the genome
            accession number. This becomes part of the sequence id printed
            in the protein FASTA output.
        @param duplicationPolicy: A C{str} indicating what to do if a
            to-be-inserted accession number is already present in the database.
            "error" results in a ValueError being raised, "ignore" means ignore
            the duplicate. It should also be possible to update (i.e., replace)
            but that is not supported yet.
        @param logfp: If not C{None}, a file pointer to write verbose
            progress output to.
        @raise DatabaseDuplicationError: If a duplicate accession number is
            encountered and C{duplicationPolicy} is 'error'.
        @return: A C{tuple} containing three C{int}s: the number of genome
            sequences examined (for potential addition), the number of genomes
            actually added, and the total number of proteins added.
        """
        assert self.SEQUENCE_ID_SEPARATOR not in proteinSource, (
            "proteinSource cannot contain %r as that is used as a separator."
            % self.SEQUENCE_ID_SEPARATOR
        )

        assert self.SEQUENCE_ID_SEPARATOR not in genomeSource, (
            "genomeSource cannot contain %r as that is used as a separator."
            % self.SEQUENCE_ID_SEPARATOR
        )

        assert not (dnaOnly and rnaOnly), "dnaOnly and rnaOnly cannot both be True."

        examinedGenomeCount = addedGenomeCount = addedProteinCount = 0

        for genome in genomes:
            examinedGenomeCount += 1
            source = getSourceInfo(genome, logfp=logfp)

            if source is None:
                # The lack of a source is logged by getSourceInfo.
                continue

            genomeLength = len(str(genome.seq))

            if logfp:
                print("\n%s: %s" % (genome.id, genome.description), file=logfp)
                print("  length = %d" % genomeLength, file=logfp)
                print("  Source:", file=logfp)
                for k, v in source.items():
                    print("    %s = %r" % (k, v), file=logfp)
                print("  Annotations:", file=logfp)
                for k, v in genome.annotations.items():
                    if k not in ("references", "comment", "structured_comment"):
                        print("    %s = %r" % (k, v), file=logfp)

            if minGenomeLength is not None and genomeLength < minGenomeLength:
                if logfp:
                    print("  Genome too short. Skipping.", file=logfp)
                continue

            if maxGenomeLength is not None and genomeLength > maxGenomeLength:
                if logfp:
                    print("  Genome too long. Skipping.", file=logfp)
                continue

            try:
                lineage = lineageFetcher(genome)
            except ValueError as e:
                print(
                    "ValueError calling lineage fetcher for %s (%s): %s"
                    % (genome.id, genome.description, e),
                    file=logfp,
                )
                lineage = taxonomyId = None
            else:
                taxonomyId = lineage[0][0]

            if allowedTaxonomicRanks:
                if lineage:
                    if not isAllowedTaxonomicRank(allowedTaxonomicRanks, lineage):
                        if logfp:
                            print(
                                "  %s (%s) has an unwanted taxonomic lineage. Skipping."
                                % (genome.id, genome.description),
                                file=logfp,
                            )
                        continue
                else:
                    # We are allowing only certain taxonomic ranks and we
                    # have no lineage information for this virus. We skip
                    # it because we cannot confirm that we want it.
                    if logfp:
                        print(
                            "  %s (%s) has no taxonomic lineage information. Skipping."
                            % (genome.id, genome.description),
                            file=logfp,
                        )
                    continue

            if dnaOnly:
                if not source["mol_type"].endswith("DNA"):
                    if logfp:
                        print(
                            "  %s (%s) is not a DNA virus (mol_type). Skipping."
                            % (genome.id, genome.description),
                            file=logfp,
                        )
                    continue
                # if lineage:
                #     print('  Lineage:', file=logfp)
                #     print(formatLineage(lineage, prefix='    '), file=logfp)
                #     if isDNAVirus(lineage):
                #         if logfp:
                #             print('  %s (%s) is a DNA virus.' %
                #                   (genome.id, genome.description),
                #                   file=logfp)
                #     else:
                #         if logfp:
                #             print('  %s (%s) is not a DNA virus.' %
                #                   (genome.id, genome.description),
                #                   file=logfp)
                #         continue
                # else:
                #     print('Could not look up taxonomy lineage for %s (%s). '
                #           'Cannot confirm as DNA.' %
                #           (genome.id, genome.description), file=logfp)
                #     continue

            if rnaOnly:
                if not source["mol_type"].endswith("RNA"):
                    if logfp:
                        print(
                            "  %s (%s) is not a RNA virus (mol_type). Skipping."
                            % (genome.id, genome.description),
                            file=logfp,
                        )
                    continue
                # if lineage:
                #     print('  Lineage:', file=logfp)
                #     print(formatLineage(lineage, prefix='    '), file=logfp)
                #     if isRNAVirus(lineage):
                #         if logfp:
                #             print('  %s (%s) is an RNA virus.' %
                #                   (genome.id, genome.description),
                #                   file=logfp)
                #     else:
                #         if logfp:
                #             print('  %s (%s) is not an RNA virus. Skipping.'
                #                   % (genome.id, genome.description),
                #                   file=logfp)
                #         continue
                # else:
                #     print('Could not look up taxonomy lineage for %s (%s). '
                #           'Cannot confirm as RNA. Skipping.' %
                #           (genome.id, genome.description), file=logfp)
                #     continue

            if excludeFungusOnlyViruses:
                if lineage is None:
                    print(
                        "  Could not look up taxonomy lineage for %s "
                        "(%s). Cannot confirm as fungus-only virus."
                        % (genome.id, genome.description),
                        file=logfp,
                    )
                else:
                    if taxonomyDatabase.isFungusOnlyVirus(lineage, genome.description):
                        if logfp:
                            print(
                                "  %s (%s) is a fungus-only virus. Skipping."
                                % (genome.id, genome.description),
                                file=logfp,
                            )
                        continue
                    else:
                        if logfp:
                            print(
                                "  %s (%s) is not a fungus-only virus."
                                % (genome.id, genome.description),
                                file=logfp,
                            )

            if excludePlantOnlyViruses:
                if lineage is None:
                    print(
                        "  Could not look up taxonomy lineage for %s "
                        "(%s). Cannot confirm as plant-only virus. "
                        "Not skipping." % (genome.id, genome.description),
                        file=logfp,
                    )
                else:
                    if taxonomyDatabase.isPlantOnlyVirus(lineage, genome.description):
                        if logfp:
                            print(
                                "  %s (%s) is a plant-only virus. Skipping."
                                % (genome.id, genome.description),
                                file=logfp,
                            )
                        continue
                    else:
                        if logfp:
                            print(
                                "  %s (%s) is not a plant-only virus."
                                % (genome.id, genome.description),
                                file=logfp,
                            )

            if excludeExclusiveHosts:
                if taxonomyId is None:
                    print(
                        "  Could not find taxonomy id for %s (%s). "
                        "Cannot exclude due to exclusive host criteria."
                        % (genome.id, genome.description),
                        file=logfp,
                    )
                else:
                    hosts = taxonomyDatabase.hosts(taxonomyId)
                    if hosts is None:
                        print(
                            "  Could not find hosts for %s (%s). Cannot "
                            "exclude due to exclusive host criteria."
                            % (genome.id, genome.description),
                            file=logfp,
                        )
                    else:
                        if len(hosts) == 1:
                            host = hosts.pop()
                            if host in excludeExclusiveHosts:
                                print(
                                    "  Skipping %s (%s) due to exclusive "
                                    "host criteria (infects only %s hosts)."
                                    % (genome.id, genome.description, host),
                                    file=logfp,
                                )
                                continue

            proteinCount = len(list(self._genomeProteins(genome)))

            if self.addGenome(
                genome,
                source,
                taxonomyId,
                proteinCount,
                databaseName,
                duplicationPolicy=duplicationPolicy,
                logfp=logfp,
            ):
                self.addProteins(
                    genome,
                    source,
                    proteinSource=proteinSource,
                    genomeSource=genomeSource,
                    duplicationPolicy=duplicationPolicy,
                    logfp=logfp,
                )

                addedProteinCount += proteinCount
                addedGenomeCount += 1

                print(
                    "  Added %s (%s) with %d protein%s to database."
                    % (
                        genome.id,
                        genome.description,
                        proteinCount,
                        "" if proteinCount == 1 else "s",
                    ),
                    file=logfp,
                )

        return examinedGenomeCount, addedGenomeCount, addedProteinCount

    def addGenome(
        self,
        genome,
        source,
        taxonomyId,
        proteinCount,
        databaseName,
        duplicationPolicy="error",
        logfp=None,
    ):
        """
        Add information about a genome to the genomes table.

        @param genome: A GenBank genome record, as parsed by SeqIO.parse
        @param source: A C{dict} containing genome source information, as
            returned by C{getSourceInfo}.
        @param taxonomyId: Either an C{int} taxonomy id or C{None} if the
            genome taxonomy could not be looked up.
        @param proteinCount: The C{int} number of proteins in the genome.
        @param databaseName: A C{str} indicating the database the records
            in C{filename} came from (e.g., 'refseq' or 'RVDB').
        @param duplicationPolicy: A C{str} indicating what to do if a
            to-be-inserted accession number is already present in the database.
            "error" results in a ValueError being raised, "ignore" means ignore
            the duplicate. It should also be possible to update (i.e., replace)
            but that is not supported yet.
        @param logfp: If not C{None}, a file pointer to write verbose
            progress output to.
        @raise DatabaseDuplicationError: If a duplicate accession number is
            encountered and C{duplicationPolicy} is 'error'.
        @return: C{True} if the genome was added, else C{False}.
        """
        sequence = str(genome.seq)

        try:
            self._connection.execute(
                "INSERT INTO genomes(accession, organism, name, sequence, "
                "length, proteinCount, host, note, taxonomyId, databaseName) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    genome.id,
                    source["organism"],
                    genome.description,
                    sequence,
                    len(sequence),
                    proteinCount,
                    source["host"],
                    source.get("note"),
                    taxonomyId,
                    databaseName,
                ),
            )
        except sqlite3.IntegrityError as e:
            if str(e).find("UNIQUE constraint failed") > -1:
                if duplicationPolicy == "error":
                    raise DatabaseDuplicationError(
                        "Genome information for %r already present in "
                        "database: %s" % (genome.id, e)
                    )
                elif duplicationPolicy == "ignore":
                    if logfp:
                        print(
                            "Genome information for %r already present in "
                            "database. Ignoring: %s" % (genome.id, e),
                            file=logfp,
                        )
                    return False
                else:
                    raise NotImplementedError(
                        "Unknown duplication policy (%s) found when "
                        "attempting to insert genome information for %s."
                        % (duplicationPolicy, genome.id)
                    )
            else:
                raise
        else:
            return True

    def addProteins(
        self,
        genome,
        source,
        proteinSource="GENBANK",
        genomeSource="GENBANK",
        duplicationPolicy="error",
        logfp=None,
    ):
        """
        Add proteins from a Genbank genome record to the proteins database and
        write out their sequences to the proteins FASTA file (in
        C{self._fastaFp}).

        @param genome: Either a GenBank genome record, as parsed by
            C{SeqIO.parse} or a C{_Genome} instance (which behaves like the
            former).
        @param source: A C{dict} containing genome source information, as
            returned by C{getSourceInfo}.
        @param proteinSource: A C{str} giving the source of the protein
            accession number. This becomes part of the sequence id printed
            in the protein FASTA output.
        @param genomeSource: A C{str} giving the source of the genome
            accession number. This becomes part of the sequence id printed
            in the protein FASTA output.
        @param duplicationPolicy: A C{str} indicating what to do if a
            to-be-inserted accession number is already present in the database.
            "error" results in a ValueError being raised, "ignore" means ignore
            the duplicate. It should also be possible to update (i.e., replace)
            but that is not supported yet.
        @param logfp: If not C{None}, a file pointer to write verbose
            progress output to.
        @raise DatabaseDuplicationError: If a duplicate accession number is
            encountered and C{duplicationPolicy} is 'error'.
        """
        genomeLen = len(genome.seq)

        for fInfo in self._genomeProteins(genome, logfp=logfp):
            # Write FASTA for the protein.
            seqId = self.SEQUENCE_ID_SEPARATOR.join(
                (
                    self.SEQUENCE_ID_PREFIX,
                    proteinSource,
                    fInfo["proteinId"],
                    genomeSource,
                    genome.id,
                    fInfo["product"],
                )
            )

            print(
                ">%s [%s]\n%s" % (seqId, source["organism"], fInfo["translation"]),
                file=self._fastaFp,
            )

            self.addProtein(
                fInfo["proteinId"],
                genome.id,
                fInfo["translation"],
                fInfo["featureLocation"],
                fInfo["forward"],
                fInfo["circular"],
                fInfo["ranges"].distinctRangeCount(genomeLen),
                gene=fInfo["gene"],
                note=fInfo["note"],
                product=fInfo["product"],
                duplicationPolicy=duplicationPolicy,
                logfp=logfp,
            )

    def addProtein(
        self,
        accession,
        genomeAccession,
        sequence,
        offsets,
        forward,
        circular,
        rangeCount,
        gene=None,
        note=None,
        product=None,
        duplicationPolicy="error",
        logfp=None,
    ):
        """
        Add information about a protein to the proteins table.

        @param accession: A C{str} protein accession id.
        @param genomeAccession: A C{str} genome accession id (the genome to
            which this protein belongs).
        @param sequence: A C{str} protein amino acid sequence.
        @param offsets: A C{str} describing the offsets of the protein in the
            genome (as obtained from C{SeqIO.parse} on a GenBank file).
        @param forward: A C{bool}, C{True} if the protein occurs on the
            forward strand of the genome, C{False} if on the complement strand.
            Note that this is converted to an C{int} in the database.
        @param circular: A C{bool}, C{True} if the protein crosses the genome
            boundary and is therefore circular, C{False} if not. Note that
            this is converted to an C{int} in the database.
        @param rangeCount: The C{int} number of ranges (regions) the protein
            comes from in the genome.
        @param gene: A C{str} gene name, or C{None} if no gene is known.
        @param note: A C{str} note about the protein, or C{None}.
        @param product: A C{str} description of the protein product (e.g.,
            "putative replication initiation protein"), or C{None}.
        @param duplicationPolicy: A C{str} indicating what to do if a
            to-be-inserted accession number is already present in the database.
            "error" results in a ValueError being raised, "ignore" means ignore
            the duplicate. It should also be possible to update (i.e., replace)
            but that is not supported yet.
        @param logfp: If not C{None}, a file pointer to write verbose
            progress output to.
        @raise DatabaseDuplicationError: If a duplicate accession number is
            encountered and C{duplicationPolicy} is 'error'.
        """
        try:
            self._connection.execute(
                "INSERT INTO proteins("
                "accession, genomeAccession, sequence, length, offsets, "
                "forward, circular, rangeCount, gene, note, product) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    accession,
                    genomeAccession,
                    sequence,
                    len(sequence),
                    offsets,
                    int(forward),
                    int(circular),
                    rangeCount,
                    gene,
                    note,
                    product,
                ),
            )
        except sqlite3.IntegrityError as e:
            if str(e).find("UNIQUE constraint failed") > -1:
                if duplicationPolicy == "error":
                    raise DatabaseDuplicationError(
                        "Protein information for %r already present in "
                        "database." % accession
                    )
                elif duplicationPolicy == "ignore":
                    if logfp:
                        print(
                            "Protein information for %r already present in "
                            "database. Ignoring." % accession,
                            file=logfp,
                        )
                else:
                    raise NotImplementedError(
                        "Unknown duplication policy (%s) found when "
                        "attempting to insert protein information for %s."
                        % (duplicationPolicy, accession)
                    )
            else:
                raise
        else:
            if logfp:
                print(
                    "    Protein %s: genome=%s product=%s"
                    % (accession, genomeAccession, product),
                    file=logfp,
                )

    def _genomeProteins(self, genome, logfp=None):
        """
        Get proteins (CDS features) that we can process from a genome, along
        with information extracted from each.

        @param genome: A GenBank genome record, as parsed by SeqIO.parse
        @param logfp: If not C{None}, a file pointer to write verbose
            progress output to.
        @return: A generator yielding feature info C{dict}s as returned by
            C{getCDSInfo}.
        """
        for feature in genome.features:
            if feature.type == "CDS":
                featureInfo = getCDSInfo(genome, feature)
                if featureInfo:
                    yield featureInfo

    def close(self):
        """
        Create indices on the accesssion ids and close the connection.
        """
        cur = self._connection.cursor()
        cur.execute(
            "CREATE UNIQUE INDEX IF NOT EXISTS protein_idx ON proteins(accession)"
        )
        cur.execute(
            "CREATE UNIQUE INDEX IF NOT EXISTS genomes_idx ON genomes(accession)"
        )
        self._connection.commit()
        self._connection.close()
        self._connection = None

    def __enter__(self):
        return self

    def __exit__(self, excType, excValue, traceback):
        self.close()


class SqliteIndex:
    """
    Provide lookup access to an Sqlite3 database holding information about
    proteins and the genomes they come from.

    @param dbFilenameOrConnection: Either a C{str} file name containing an
        sqlite3 database as created by C{SqliteIndexWriter} or an already
        open connection to such a database. Note that an already open
        connection will not be closed by self.close().
    @param lookupCacheSize: The C{int} size of the memoization cache
        for the protein and genome lookup functions (each has its own
        memoization cache).
    """

    PROTEIN_ACCESSION_FIELD = 2
    GENOME_ACCESSION_FIELD = 4

    def __init__(self, dbFilenameOrConnection, lookupCacheSize=1024):
        if isinstance(dbFilenameOrConnection, str):
            self._connection = sqliteConnect(dbFilenameOrConnection)
            self._closeConnection = True
        else:
            self._connection = dbFilenameOrConnection
            self._closeConnection = False
        self._connection.row_factory = sqlite3.Row
        self._proteinCache = LRUCache(maxsize=lookupCacheSize)
        self._genomeCache = LRUCache(maxsize=lookupCacheSize)

    def genomeAccession(self, id_):
        """
        Get the genome accession info from a sequence id.

        @param id_: A C{str} sequence id in the form
            'civ|GENBANK|%s|GENBANK|%s|%s [%s]' where the genome accession
            is in the fifth '|'-separated field.
        @raise IndexError: If C{id_} does not have enough |-separated fields.
        @return: The C{str} accession number.
        """
        return id_.split("|", self.GENOME_ACCESSION_FIELD + 1)[
            self.GENOME_ACCESSION_FIELD
        ]

    def proteinAccession(self, id_):
        """
        Get the protein accession info from a sequence id.

        @param id_: A C{str} sequence id in the form
            'civ|GENBANK|%s|GENBANK|%s|%s [%s]' where the protein accession
            is in the third '|'-separated field.
        @raise IndexError: If C{id_} does not have enough |-separated fields.
        @return: The C{str} accession number.
        """
        return id_.split("|", self.PROTEIN_ACCESSION_FIELD + 1)[
            self.PROTEIN_ACCESSION_FIELD
        ]

    @cachedmethod(attrgetter("_genomeCache"))
    def _findGenome(self, accession):
        """
        Find info about a genome, given an accession number.

        @param accession: A C{str} accession number.
        @return: A C{dict} with keys corresponding to the names of the columns
            in the genomes database table, else C{None} if C{id_} cannot be
            found.
        """
        cur = self.execute("SELECT * FROM genomes WHERE accession = ?", (accession,))
        row = cur.fetchone()

        if row:
            result = dict(row)
            # TODO: the following line can be removed, I think.
            result["accession"] = accession
            return result

    def findGenome(self, id_):
        """
        Find info about a genome, given a sequence id.

        @param id_: A C{str} sequence id. This is either of the form
            'civ|GENBANK|%s|GENBANK|%s|%s [%s]' where the genome id is in the
            5th '|'-delimited field, or else is the nucleotide sequence
            accession number as already extracted.
        @return: A C{dict} with keys corresponding to the names of the columns
            in the genomes database table, else C{None} if C{id_} cannot be
            found.
        """
        try:
            accession = self.genomeAccession(id_)
        except IndexError:
            accession = id_

        return self._findGenome(accession)

    @cachedmethod(attrgetter("_proteinCache"))
    def _findProtein(self, accession):
        """
        Find info about a protein, given an accession number.

        @param accession: A C{str} accession number.
        @return: A C{dict} with keys corresponding to the names of the columns
            in the proteins database table, else C{None} if C{id_} cannot be
            found.
        """
        cur = self.execute("SELECT * FROM proteins WHERE accession = ?", (accession,))
        row = cur.fetchone()
        if row:
            result = dict(row)
            result["forward"] = bool(result["forward"])
            result["circular"] = bool(result["circular"])
            result["length"] = int(result["length"])
            # TODO: the following line can be removed, I think.
            result["accession"] = accession
            return result

    def findProtein(self, id_):
        """
        Find info about a protein, given a sequence id.

        @param id_: A C{str} sequence id. This is either of the form
            'civ|GENBANK|%s|GENBANK|%s|%s [%s]' where the protein id is in the
            3rd '|'-delimited field, or else is the protein accession number as
            already extracted.
        @return: A C{dict} with keys corresponding to the names of the columns
            in the proteins database table, else C{None} if C{id_} cannot be
            found.
        """
        try:
            accession = self.proteinAccession(id_)
        except IndexError:
            accession = id_

        return self._findProtein(accession)

    def _yieldProteins(self, rows, cur):
        """
        Helper function for self.findProteinsForGenome.

        @param rows: A C{list} of protein database lookup results.
        @param cur: An sqlite3 cursor.
        @return: A generator that yields C{dict}s with keys corresponding to
            the names of the columns in the proteins database table.
        """
        while rows:
            for row in rows:
                result = dict(row)
                result["forward"] = bool(result["forward"])
                result["circular"] = bool(result["circular"])
                result["length"] = int(result["length"])
                yield result
            rows = cur.fetchmany()

    def findProteinsForGenome(self, id_):
        """
        Find all proteins for a genome id.

        @param id_: A C{str} sequence id. This is either of the form
            'civ|GENBANK|%s|GENBANK|%s|%s [%s]' where the genome id is in the
            5th '|'-delimited field, or else is the nucleotide sequence
            accession number as already extracted.
        @return: A generator that yields C{dict}s with keys corresponding to
            the names of the columns in the proteins database table, else
            C{None} if C{id_} cannot be found.
        """
        try:
            accession = self.genomeAccession(id_)
        except IndexError:
            accession = id_

        cur = self.execute(
            "SELECT * FROM proteins WHERE genomeAccession = ?", (accession,)
        )

        rows = cur.fetchmany()
        if rows:
            return self._yieldProteins(rows, cur)

    def execute(self, query, *args):
        """
        Execute an SQL statement. See
        https://docs.python.org/3.5/library/sqlite3.html#sqlite3.Cursor.execute
        for full argument details.

        @param query: A C{str} SQL query.
        @param args: Additional arguments (if any) to pass to the sqlite3
            execute command.
        @return: An sqlite3 cursor.
        """
        cur = self._connection.cursor()
        cur.execute(query, *args)
        return cur

    def proteinCount(self):
        """
        How many proteins are in the database?

        @return: An C{int} protein count.
        """
        cur = self.execute("SELECT COUNT(1) FROM proteins")
        return int(cur.fetchone()[0])

    def genomeCount(self):
        """
        How many genomes are in the database?

        @return: An C{int} genome count.
        """
        cur = self.execute("SELECT COUNT(1) FROM genomes")
        return int(cur.fetchone()[0])

    def close(self):
        """
        Close the database connection (if we opened it).
        """
        if self._closeConnection:
            self._connection.close()
        self._connection = None

    def __enter__(self):
        return self

    def __exit__(self, excType, excValue, traceback):
        self.close()
