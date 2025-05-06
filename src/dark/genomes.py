from collections import Counter

from dark.errors import NoSuchGenomeError
from dark.genbank import GenomeRanges
from dark.reads import DNARead
from dark.sam import samfile


class GenomeProteinInfo:
    """
    Hold information about the proteins in a genome and how they are matched
    by reads in SAM files.

    @param accession: The C{str} accession number of a genome.
    @param proteinGenomeDB: A L{dark.civ.proteins.SqliteIndex} instance.
    @param checkTranslations: If C{True}, check that the protein sequences
        that area supposed to have come from the genome can be obtained by
        translating the corresponding region of the genome.
    """

    def __init__(self, genomeAccession, proteinGenomeDB, checkTranslations=True):
        self.genomeAccession = genomeAccession
        self.proteinGenomeDB = proteinGenomeDB
        # self.proteins is keyed by protein accession number.
        self.proteins = {}
        self.coveredProteins = set()
        # self.offsets is keyed by genome offset, values are dicts that
        # contain a list of protein accession numbers that overlap that
        # offset and a set of read ids (if any) that match at that offset.
        # The offsets keys are only those that correspond to one or more
        # proteins in the genome.
        self.offsets = {}
        # self.coveredOffsetCount holds the read counts for all offsets covered
        # by reads, regardless of whether the offsets correspond to proteins or
        # not.
        self.coveredOffsetCount = Counter()
        self.samFiles = []
        self.readIdsMatchingGenome = set()

        self.genome = proteinGenomeDB.findGenome(genomeAccession)
        if self.genome is None:
            raise NoSuchGenomeError(
                "Reference %r not found in protein/genome "
                "database." % genomeAccession
            )

        for protein in proteinGenomeDB.findProteinsForGenome(genomeAccession):
            proteinAccession = protein["accession"]
            self.proteins[proteinAccession] = protein

            ranges = GenomeRanges(protein["offsets"]).ranges
            # print('Protein accession', proteinAccession)
            # print(ranges)

            for start, stop, forward in ranges:
                for offset in range(start, stop):
                    if offset not in self.offsets:
                        self.offsets[offset] = {
                            "proteinAccessions": set(),
                            "readIds": set(),
                        }
                    self.offsets[offset]["proteinAccessions"].add(proteinAccession)

            if checkTranslations:
                self._checkTranslation(self.genome, ranges, protein)

    def _checkTranslation(self, genome, ranges, protein):
        """
        Make sure all protein sequences supposed to be in the genome can in
        fact be obtained by translating the genome.

        @param genome: A C{dict} with genome information from our sqlite3
            protein/genome database, as returned by
            C{dark.civ.proteins.SqliteIndex.findGenome.
        @param ranges: A C{list} of (start, stop, forward) nucleotide ranges
            for the protein in the genome.
        @param protein: A C{dict} with protein information from our sqlite3
            protein/genome database, as returned by
            C{dark.civ.proteins.SqliteIndex.findProtein.
        """
        proteinSequence = protein["sequence"] + "*"

        # print('protein name', protein['product'], 'ranges', ranges)
        sequence = "".join(
            [genome["sequence"][start:stop] for (start, stop, _) in ranges]
        )

        genomeRead = DNARead("id", sequence)
        translations = list(genomeRead.translations())
        index = 0 if protein["forward"] else 3
        if translations[index].sequence != proteinSequence:
            # TODO: improve this error to show what actually went wrong.
            raise ValueError("Could not translate genome range to get protein sequence")

    def addSAM(self, filename, filterAlignment=None):
        """
        Read a SAM file and add information about the reads that match our
        reference id.

        @param filename: A C{str} SAM filename.
        @param filterAlignment: A 1-argument function to be used for filtering
            reads in the SAM file. If C{None}, all alignments will be examined.
        """
        self.samFiles.append(filename)
        referenceId = self.genomeAccession
        with samfile(filename) as sam:
            for column in sam.pileup():
                for read in column.pileups:
                    alignment = read.alignment
                    if alignment.reference_name == referenceId and (
                        filterAlignment is None or filterAlignment(alignment)
                    ):
                        readId = alignment.query_name
                        self.readIdsMatchingGenome.add(readId)
                        offset = column.reference_pos
                        self.coveredOffsetCount[offset] += 1

                        try:
                            offsetInfo = self.offsets[offset]
                        except KeyError:
                            pass
                        else:
                            # This offset corresponds to one or more proteins.
                            self.coveredProteins.update(offsetInfo["proteinAccessions"])
                            offsetInfo["readIds"].add(readId)

    def proteinCoverageInfo(self, proteinAccession, minReadOffsetCount=None):
        """
        Calculate coverage information for a protein.

        @param proteinAccession: A C{str} accession number.
        @param minReadOffsetCount: An C{int}, specifying the minimum number of
            reads offsets that must overlap the protein for the read to be
            considered as sufficiently intersecting the protein. Use this to
            prevent reads that just overlap the protein in a very small number
            offsets from being counted. Or C{None} to indicate that no such
            filtering should be applied.
        @raises KeyError: If C{proteinAccession} is not known.
        @return: A C{dict} containing
                * the number of covered offsets,
                * the total number of read bases that cover the protein,
                * the protein length (in nucleotides), and
                * the set of all matching read ids.
            See below for the dictionary keys.
        """
        protein = self.proteins[proteinAccession]
        coveredOffsets = 0
        totalBases = 0
        allReadIds = set()
        offsetsSeen = set()
        proteinLength = 0

        if minReadOffsetCount is not None and minReadOffsetCount < 2:
            # A minimum of zero or one is equivalent to not giving a value.
            minReadOffsetCount = None

        if minReadOffsetCount:
            readOffsetCounts = Counter()

        proteinRanges = GenomeRanges(protein["offsets"]).ranges

        # Do an initial pass across all the offsets of the protein to see
        # which reads intersect and where. We will then do a second pass in
        # which we ignore reads that do not sufficiently overlap.
        for start, stop, forward in proteinRanges:
            proteinLength += stop - start
            for offset in range(start, stop):
                assert offset not in offsetsSeen
                offsetsSeen.add(offset)
                readIds = self.offsets[offset]["readIds"]
                if readIds and minReadOffsetCount:
                    readOffsetCounts.update(readIds)

        # Sanity check that the sum of the range lengths is the same as the
        # overall length given in the database.
        #
        # The +3 in the following is because the database holds the AA
        # length, not including the stop codon. But the database range
        # covers the stop codon.
        dbProteinLength = self.proteins[proteinAccession]["length"] * 3 + 3
        if proteinLength != dbProteinLength:
            raise ValueError(
                "Sum of protein database ranges (%d) does not agree with "
                "database protein length (%d) for protein %s!"
                % (proteinLength, dbProteinLength, proteinAccession)
            )

        # If we are not reporting reads whose overlapping offset count is
        # too low, make a set of such reads.
        if minReadOffsetCount:
            unwanted = set(
                readId
                for readId in readOffsetCounts
                if readOffsetCounts[readId] < minReadOffsetCount
            )
        else:
            unwanted = set()

        # Second pass, in which we ignore unwanted (i.e., insufficiently
        # overlapping) reads.
        for start, stop, forward in proteinRanges:
            for offset in range(start, stop):
                readIds = set(self.offsets[offset]["readIds"]) - unwanted
                if readIds:
                    allReadIds.update(readIds)
                    coveredOffsets += 1
                    totalBases += len(readIds)

        return {
            "coveredOffsets": coveredOffsets,
            "totalBases": totalBases,
            "ntLength": proteinLength,
            "readIds": allReadIds,
        }

    def readIdsForAllProteins(self):
        """
        Get the set of read ids for reads matching any protein at all.

        @return: A C{set} of C{str} read ids.
        """
        result = set()
        for offsetInfo in self.offsets.values():
            result.update(offsetInfo["readIds"])
        return result
