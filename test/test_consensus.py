from unittest import TestCase, skipUnless

from dark.consensus import consensusFromBAM, ConsensusError
from dark.reads import DNARead
from dark.sam import (
    samtoolsInstalled,
    UnequalReferenceLengthError,
    UnknownReference,
    UnspecifiedReference,
)

from .bam import makeBAM, REF_ID


@skipUnless(samtoolsInstalled(), "samtools is not installed")
class TestReferenceErrors(TestCase):
    """
    Test expected errors.
    """

    def testUnknownStrategy(self):
        """
        If an unknown strategy is passed, ValueError must be raised.
        """
        template = ("ACGTTCCG",)

        with makeBAM(template) as (fastaFilename, bamFilename):
            strategy = "xxx"
            error = rf"^Unknown consensus strategy {strategy!r}\.$"
            self.assertRaisesRegex(
                ConsensusError,
                error,
                consensusFromBAM,
                bamFilename,
                strategy=strategy,
                quiet=True,
                referenceFasta=fastaFilename,
            )

    def testUnknownReferenceId(self):
        """
        If a reference id that is not in the BAM file is passed,
        UnknownReference must be raised.
        """
        template = ("ACGTTCCG",)

        with makeBAM(template) as (fastaFilename, bamFilename):
            badId = "bad-id"
            error = (
                rf"^BAM file {str(bamFilename)!r} does not mention a "
                rf"reference with id {badId!r}\. Known references are: "
                rf"ref-id\.$"
            )
            self.assertRaisesRegex(
                UnknownReference,
                error,
                consensusFromBAM,
                bamFilename,
                bamId=badId,
                referenceFasta=fastaFilename,
            )

    def testUnknownReference(self):
        """
        If a reference sequence is passed and its id is not in the BAM
        file, UnknownReference must be raised.
        """
        template = ("ACGTTCCG",)

        UNKNOWN = "unknown-id"
        with makeBAM(template) as (fastaFilename, bamFilename):
            error = (
                rf"^BAM file {str(bamFilename)!r} does not mention a "
                rf"reference with id {UNKNOWN!r}\. Known "
                rf"references are: ref-id\.$"
            )
            self.assertRaisesRegex(
                UnknownReference,
                error,
                consensusFromBAM,
                bamFilename,
                bamId=UNKNOWN,
                referenceFasta=fastaFilename,
            )

    def testUnknownFastaReference(self):
        """
        If a reference sequence id is passed and it is not in the FASTA file,
        UnknownReference must be raised.
        """
        template = ("ACGTTCCG",)

        UNKNOWN = "unknown-id"
        with makeBAM(template) as (fastaFilename, bamFilename):
            error = (
                rf"^No sequence with id {UNKNOWN!r} found in "
                rf"{str(fastaFilename)!r}\.$"
            )
            self.assertRaisesRegex(
                UnknownReference,
                error,
                consensusFromBAM,
                bamFilename,
                fastaId=UNKNOWN,
                referenceFasta=fastaFilename,
            )

    def testUnspecifiedReference(self):
        """
        If a reference is not specified and the BAM file mentions more than
        one, UnspecifiedReference must be raised.
        """
        template = ("ACGTTCCG",)

        bamReferences = [DNARead("ref-1", template[0]), DNARead("ref-2", "AA")]
        fastaReferences = [DNARead("ref-3", "AAA")]
        with makeBAM(
            template, bamReferences=bamReferences, fastaReferences=fastaReferences
        ) as (fastaFilename, bamFilename):
            error = (
                r"^Could not infer a BAM reference. Available references are: "
                r"ref-1, ref-2\.$"
            )
            self.assertRaisesRegex(
                UnspecifiedReference,
                error,
                consensusFromBAM,
                bamFilename,
                referenceFasta=fastaFilename,
            )

    def testReferenceOfWrongLength(self):
        """
        If a reference sequence is passed and it does not have the right
        length, UnequalReferenceLengthError must be raised.
        """
        template = ("ACGTTCCG",)

        fastaReferences = [DNARead(REF_ID, "AA")]
        with makeBAM(template, fastaReferences=fastaReferences) as (
            fastaFilename,
            bamFilename,
        ):
            error = (
                rf"^Reference FASTA sequence {REF_ID!r} has length 2, but the "
                rf"BAM reference {REF_ID!r} has length {len(template[0])}\.$"
            )
            self.assertRaisesRegex(
                UnequalReferenceLengthError,
                error,
                consensusFromBAM,
                bamFilename,
                referenceFasta=fastaFilename,
                quiet=True,
            )

    def testEmptyReferenceFile(self):
        """
        If a reference file is passed but is empty, UnknownReference must be
        raised.
        """
        template = ("ACGTTCCG",)

        fastaReferences = []
        with makeBAM(template, fastaReferences=fastaReferences) as (
            fastaFilename,
            bamFilename,
        ):
            error = (
                rf"^The FASTA reference file {str(fastaFilename)!r} contained "
                rf"no sequences\.$"
            )
            self.assertRaisesRegex(
                UnknownReference,
                error,
                consensusFromBAM,
                bamFilename,
                referenceFasta=fastaFilename,
            )


class _Mixin:
    """
    Common (i.e., with and without considering FASTQ quality values) tests for
    consensuses making.
    """

    def testNoReadsReference(self):
        """
        If no reads are present and resolution of no-coverage bases is the
        reference sequence, the reference should be returned.
        """
        template = ("ACGTTCCG",)

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                template[0],
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testNoReadsReferenceFromId(self):
        """
        The reference can be passed using just its id (as opposed to a full
        reference Read instance). A string of 'N's is returned because there
        are no reads and we pass noCoverage='N'.
        """
        template = ("ACGTTCCG",)

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "N" * len(template[0]),
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    bamId="ref-id",
                    noCoverage="N",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testReferenceAndReferenceIdNotGiven(self):
        """
        If the reference to use is not given, the single reference in the BAM
        file should be used. A string of 'N's is returned because there
        are no reads and we pass noCoverage='N'.
        """
        template = ("ACGTTCCG",)

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "N" * len(template[0]),
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    noCoverage="N",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testNoReadsN(self):
        """
        If no reads are present and resolution of no-coverage bases is 'N'
        (or '?', etc) a sequence of Ns (or ?s, etc) should be returned.
        """
        template = ("ACGTTCCG",)

        with makeBAM(template) as (fastaFilename, bamFilename):
            for char in "N?":
                self.assertEqual(
                    char * len(template[0]),
                    consensusFromBAM(
                        bamFilename,
                        quiet=True,
                        referenceFasta=fastaFilename,
                        noCoverage=char,
                        ignoreQuality=self.ignoreQuality,
                    ).sequence,
                )

    def testLowReadsReference(self):
        """
        If fewer reads than needed are present and resolution of low-coverage
        sites is the reference sequence, the reference should be returned.
        """
        template = (
            "ACGTTCCG",
            "  GTT",
            "  ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                template[0],
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    minCoverage=2,
                    lowCoverage="reference",
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testLowReadsCharNoCoverageConsensus(self):
        """
        If fewer reads than needed are present and resolution of low-coverage
        sites is 'N' (or '? etc), then Ns (or ?s, etc) should be returned in
        the low coverage sites, and the reference in the sites with no
        coverage.
        """
        template = (
            "ACGTTCCG",
            "  GTT",
            "  ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            for char in "N?":
                self.assertEqual(
                    "AC" + char * 3 + "CCG",
                    consensusFromBAM(
                        bamFilename,
                        quiet=True,
                        referenceFasta=fastaFilename,
                        minCoverage=2,
                        lowCoverage=char,
                        noCoverage="reference",
                        ignoreQuality=self.ignoreQuality,
                    ).sequence,
                )

    def testLowReadsCharNoCoverageX(self):
        """
        If fewer reads than needed are present and resolution of low-coverage
        sites is 'N' (or '? etc), then Ns (or ?s, etc) should be returned in
        the low coverage sites, and (for example) 'X' in the sites with no
        coverage.
        """
        template = (
            "ACGTTCCG",
            "  GTT",
            "  ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            for char in "N?":
                self.assertEqual(
                    "XX" + char * 3 + "XXX",
                    consensusFromBAM(
                        bamFilename,
                        quiet=True,
                        referenceFasta=fastaFilename,
                        noCoverage="X",
                        lowCoverage=char,
                        minCoverage=2,
                        ignoreQuality=self.ignoreQuality,
                    ).sequence,
                )

    def testOneReadMatchingPartOfTheReference(self):
        """
        If one read is present and it matches part of the reference and
        resolution of no-coverage bases is the reference sequence, the
        reference should be returned.
        """
        template = (
            "ACGTTCCG",
            "  GTT",
            "  ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                template[0],
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOneReadDifferingFromPartOfTheReference(self):
        """
        If one read is present and it differs from part of the reference and
        resolution of no-coverage bases is the reference sequence, the
        expected hybrid of the read and the reference should be returned.
        """
        template = (
            "ACGTTCCG",
            "  AAA",
            "  ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACAAACCG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoReadsDifferingFromPartOfTheReferenceSomeLowCoverage(self):
        """
        If two reads are present and they differ from part of the reference and
        resolution of no-coverage bases is the reference sequence, the
        expected hybrid of the read and the reference should be returned, with
        the part of the reads that has insufficient coverage returning the
        low-coverage symbol (here '+').
        """
        template = (
            "ACGTTCCG",
            "  AAA",
            "  ???",
            "  AAAT",
            "  ????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACAAA+CG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    minCoverage=2,
                    noCoverage="reference",
                    lowCoverage="+",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoReadsDifferingFromPartOfTheReferenceLowAndNoCoverage(self):
        """
        If two reads are present and they differ from part of the reference and
        resolution of no-coverage bases is 'N', the expected hybrid of the read
        and the reference should be returned, with the part of the reads that
        has insufficient coverage returning the low-coverage symbol (here '+').
        """
        template = (
            "ACGTTCCG",
            "  AAA",
            "  ???",
            "  AAAT",
            "  ????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "NNAAA+NN",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    minCoverage=2,
                    noCoverage="N",
                    lowCoverage="+",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testSimpleMajority(self):
        """
        If three reads result in a majority base at a site, that base should
        be in the consensus.
        """
        template = (
            "ACGT",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  C",
            "  ?",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACAT",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.5,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testSimpleMajorityBelowThreshold(self):
        """
        If conflicting reads at a site do not give a simple (above threshold)
        majority, the ambiguous code should be in the consensus.
        """
        template = (
            "ACGT",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  C",
            "  ?",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACMT",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testReadHasEarlierSites(self):
        """
        If a read has sites that come before the reference, they must
        be included in the consensus.
        """
        template = (
            " ACGT",
            "AA",
            "??",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "AACGT",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    includeSoftClipped=True,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testReadHasEarlierSitesNoSoftClipped(self):
        """
        If a read has sites that come before the reference, they must
        not appear in the consensus if soft-clipped bases are not included.
        """
        template = (
            " ACGT",
            "AA",
            "??",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACGT",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testReadHasLaterSites(self):
        """
        If a read has sites that come after the reference, they must
        be included in the consensus.
        """
        template = (
            "ACGT",
            "   TAA",
            "   ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACGTAA",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    includeSoftClipped=True,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testReadHasLaterSitesNoSoftClipped(self):
        """
        If a read has sites that come after the reference, they must
        not appear in the consensus if soft-clipped bases are not included.
        """
        template = (
            "ACGT",
            "   TAA",
            "   ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACGT",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testReadHasEarlierAndLaterSites(self):
        """
        If a read has sites that come before and after the reference, the sites
        must be included in the consensus.
        """
        template = (
            "  ACGT",
            "TCACGTGA",
            "????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TCACGTGA",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    includeSoftClipped=True,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testReadHasEarlierAndLaterSitesNoSoftClipped(self):
        """
        If a read has sites that come before and after the reference, the sites
        must not appear in the consensus if soft-clipped bases are not
        included.
        """
        template = (
            "  ACGT",
            "TCACGTGA",
            "????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACGT",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testReadsHaveEarlierAndLaterSites(self):
        """
        If two reads have sites that come before and after the reference, the
        sites must be included in the consensus. This is the same as the
        testReadHasEarlierAndLaterSites test above, but using two reads
        (instead of one) to give the before and after sites.
        """
        template = (
            "  ACGT",
            "TCA",
            "???",
            "     TGA",
            "     ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TCACGTGA",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    includeSoftClipped=True,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testReadsHaveEarlierAndLaterSitesNoSoftClipped(self):
        """
        If two reads have sites that come before and after the reference, the
        sites must not appear in the consensus if soft-clipped bases are not
        included. This is the same as the
        testReadHasEarlierAndLaterSitesNoSoftClipped test above, but using two
        reads (instead of one) to give the before and after sites.
        """
        template = (
            "  ACGT",
            "TCA",
            "???",
            "     TGA",
            "     ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACGT",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOneDeletionFromReference(self):
        """
        A deletion from the reference must be handled correctly.
        """
        template = (
            "CACGTG",
            " A-A",
            " ?-?",
            " A-A",
            " ?-?",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAxATG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    deletionSymbol="x",
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOneLowFrequencyDeletion(self):
        """
        A deletion from the reference that does not meet the required deletion
        frequency should not appear in the consensus.
        """
        template = (
            "CACGTG",
            " A-A",
            " ?-?",
            " AGA",
            " ???",
            " AGA",
            " ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAGATG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    deletionSymbol="x",
                    deletionThreshold=0.5,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOneDeletionFromReferenceUnmarked(self):
        """
        A deletion from the reference must be handled correctly.
        """
        template = (
            "CACGTG",
            " A-A",
            " ?-?",
            " A-A",
            " ?-?",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAATG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    deletionSymbol="",
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoDeletionsFromReference(self):
        """
        Two deletions from the reference must be handled correctly.
        """
        template = (
            "CACGTG",
            " A-A-G",
            " ?-?-?",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAxAxG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    deletionSymbol="x",
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoDeletionsFromReferenceUnmarked(self):
        """
        Two deletions from the reference must be handled correctly.
        """
        template = (
            "CACGTG",
            " A-A-G",
            " ?-?-?",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAAG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    deletionSymbol="",
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOneInsertionInReference(self):
        """
        An insertion in the reference must be handled correctly.
        """
        template = (
            "CA-CGTG",
            " AAC",
            " ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAACGTG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOneInsertionInReferenceLowFrequency(self):
        """
        An insertion in the reference must not be included if it does not
        meet the insertionCountThreshold.
        """
        template = (
            "CA-CGTG",
            " ATC",
            " ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CACGTG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    insertionCountThreshold=2,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoInsertionsInReference(self):
        """
        Two insertions in the reference must be handled correctly.
        """
        template = (
            "CA-CGT-G",
            " AAC",
            " ???",
            "     TAG",
            "     ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAACGTAG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoInsertionsAndOneUnmarkedDeletionInReference(self):
        """
        Two insertions and one deletion in the reference must be handled
        correctly.
        """
        template = (
            "CA-CGT-G",
            " AA-G",
            " ??-?",
            "     TAG",
            "     ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAAGTAG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    deletionSymbol="",
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoInsertionsAndOneDeletionInReference(self):
        """
        Two insertions and one deletion in the reference must be handled
        correctly, with the deletion marked with 'x' as requested.
        """
        template = (
            "CA-CGT-G",
            " AA-G",
            " ??-?",
            "     TAG",
            "     ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAAxGTAG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    deletionSymbol="x",
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoInsertionsAndTwoDeletionsInReference(self):
        """
        Two insertions and two deletions in the reference must be handled
        correctly and the deletions must appear in the consensus as requested.
        """
        template = (
            "CA-CGT-G",
            " AA--",
            " ??--",
            "     TAG",
            "     ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAAxxTAG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    deletionSymbol="x",
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoInsertionsAndTwoUnmarkedDeletionsInReference(self):
        """
        Two insertions and two deletions in the reference must be handled
        correctly and the deletions must not appear in the consensus if the
        deletion symbol is the empty string.
        """
        template = (
            "CA-CGT-G",
            " AA--",
            " ??--",
            "     TAG",
            "     ???",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "CAATAG",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.7,
                    deletionSymbol="",
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testGeneiousExamplesNoTie(self):
        """
        Test the no-tied counts example from
        https://assets.geneious.com/manual/2020.1/static/GeneiousManualse43.html
        """
        template = (
            "ACGT",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  G",
            "  ?",
            "  G",
            "  ?",
            "  G",
            "  ?",
            "  T",
            "  ?",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            for expected, threshold in ("A", 0.4), ("R", 0.7), ("D", 0.95):
                self.assertEqual(
                    f"AC{expected}T",
                    consensusFromBAM(
                        bamFilename,
                        quiet=True,
                        referenceFasta=fastaFilename,
                        threshold=threshold,
                        noCoverage="reference",
                        ignoreQuality=self.ignoreQuality,
                    ).sequence,
                )

    def testGeneiousExamplesTie(self):
        """
        Test the tied counts example from
        https://assets.geneious.com/manual/2020.1/static/GeneiousManualse43.html
        """
        template = (
            "ACGT",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  A",
            "  ?",
            "  G",
            "  ?",
            "  G",
            "  ?",
            "  T",
            "  ?",
            "  T",
            "  ?",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            for expected, threshold in ("A", 0.4), ("D", 0.7), ("D", 0.95):
                self.assertEqual(
                    f"AC{expected}T",
                    consensusFromBAM(
                        bamFilename,
                        quiet=True,
                        referenceFasta=fastaFilename,
                        threshold=threshold,
                        noCoverage="reference",
                        ignoreQuality=self.ignoreQuality,
                    ).sequence,
                )

    def testTwoAgreeingSoftClipsNothingBefore(self):
        """
        Test that two soft-clipped regions that agree with each other are
        in the expected result when soft-clipped bases are included.
        """
        template = (
            "        TGATCTCC",
            "AGCCAGAATGATCTCC",
            "????????????????",
            "    AGAATGATCTCC",
            "    ????????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "AGCCAGAATGATCTCC",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    includeSoftClipped=True,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoAgreeingSoftClipsNothingBeforeNoSoftClipped(self):
        """
        Test that two soft-clipped regions that agree with each other are
        ignored when soft-clipped bases are not including.
        """
        template = (
            "        TGATCTCC",
            "AGCCAGAATGATCTCC",
            "????????????????",
            "    AGAATGATCTCC",
            "    ????????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TGATCTCC",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoAgreeingInsertionsOneMatchingSomethingBefore(self):
        """
        Test that two insertions that agree with each other, with one also
        matching reference bases before the insertion give the expected result.
        """
        template = (
            "TT--------TGATCTCC",
            "TTAGCCAGAATGATCTCC",
            "??????????????????",
            "      AGAATGATCTCC",
            "      ????????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TTAGCCAGAATGATCTCC",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testTwoInsertions(self):
        """
        Test that two insertions that agree with each other, with one also
        matching reference bases before the insertion give the expected result.
        """
        template = (
            "TT--------TGAT--CTCC",
            "TTAGCCAGAATGA",
            "?????????????",
            "           GATGGCTCC",
            "           ?????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TTAGCCAGAATGATGGCTCC",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOmicronEPE214Insertion(self):
        """
        Test that an amino acid EPE sequence (here 'GAGCCAGAA') insertion
        into the SARS-CoV-2 spike nucleotide sequence at location 642 (amino
        acid location 214) works as expected.

        The nucleotide sequence below can be obtained via:

        $ ncbi-fetch-id.py MN908947.3 > MN908947.3.fasta
        $ describe-genome.py --feature S --printNtSeq < MN908947.3.fasta | \
              filter-fasta.py --keepSites 630-673 --quiet | tail -n 1

        I then inserted the 9 nucleotide sequence GAGCCAGAA before the TGAT...
        starting at position 642.
        """
        template = (
            "TAATTTAGTGCG---------TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC",
            "      AGTGCGGAGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTT",
            "      ??????????????????????????????????????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TAATTTAGTGCGGAGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOmicronEPE214InsertionRightSideExact(self):
        """
        Test that an amino acid EPE sequence (here 'GAGCCAGAA') insertion
        into the SARS-CoV-2 spike nucleotide sequence at location 642 (amino
        acid location 214) works as expected when the insertion is the very
        last part of the read, in which case it will be marked as a
        soft-clipped region.
        """
        template = (
            "TAATTTAGTGCG---------TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC",
            "      AGTGCGGAGCCAGAA",
            "      ???????????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TAATTTAGTGCGGAGCCAGAATCAGGGTTTTTCGGCTTTAGAAC",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    insertionCountThreshold=1,
                    includeSoftClipped=True,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOmicronEPE214InsertionLeftSideExactSoftClipped(self):
        """
        Test that an amino acid EPE sequence (here 'GAGCCAGAA') insertion
        into the SARS-CoV-2 spike nucleotide sequence at location 642 (amino
        acid location 214) works as expected when the insertion sequence is
        the very beginning of the read, in which case it will be marked as a
        soft-clipped region.
        """
        template = (
            "TAATTTAGTGCG---------TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC",
            "            GAGCCAGAATGATCT",
            "            ???????????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TAAGAGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    insertionCountThreshold=1,
                    includeSoftClipped=True,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOmicronEPE214PartialInsertionInTwoReadsSoftClipped(self):
        """
        The amino acid EPE sequence (here 'GAGCCAGAA') insertion
        into the SARS-CoV-2 spike nucleotide sequence at location 642 (amino
        acid location 214) must work as expected, when the reads only
        partially cover it and are therefore both soft-clipped in that region.
        They overlap in a non-agreeing way and ambiguous nucleotide codes
        result.
        """
        template = (
            "TAATTTAGTGCG---------TGATCTCCCTCA",
            "      AGTGCGGAGCCA",
            "      ????????????",
            "             AGCCAGAATGATCTCCCTCA",
            "             ????????????????????",
        )

        # This example is a bit complicated. The reference and reads will
        # align as follows (with the '-' from the reference removed to make
        # it easier to see what's going on). The final line of hyphens shows
        # the ambiguous region.
        #
        # TAATTTAGTGCGTGATCTCCCTCA
        #       AGTGCGGAGCCA
        #     AGCCAGAATGATCTCCCTCA
        #       ------------

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TAATAGMSWGMRKRRYCWCCCTCA",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    insertionCountThreshold=1,
                    includeSoftClipped=True,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOmicronEPE214PartialInsertionInTwoReads(self):
        """
        This is the same as the immediately above test
        (testOmicronEPE214PartialInsertionInTwoReadsSoftClipped), but the
        reads are now fully aligned wtih the reference so there is no soft
        clipping. The alignment of the reads is the same and so is the result.
        """
        template = (
            "TAATTTAGTGCGTGATCTCCCTCA",
            "      AGTGCGGAGCCA",
            "      ????????????",
            "    AGCCAGAATGATCTCCCTCA",
            "    ????????????????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TAATAGMSWGMRKRRYCWCCCTCA",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    insertionCountThreshold=1,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )

    def testOmicronEPE214PartialNoSoftClipped(self):
        """
        Test that a trailing part of the amino acid EPE sequence (here
        'AGCCAGAA') insertion into the SARS-CoV-2 spike nucleotide sequence
        that usually is found at location 642 (amino acid location 214) is
        excluded when it appears as a partial soft-clipped region and
        soft-clipped bases are not included.
        """
        template = (
            "        TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC",
            "AGCCAGAATGATCTCCCTCAGGGTTTTTCGGCTTT",
            "???????????????????????????????????",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "TGATCTCCCTCAGGGTTTTTCGGCTTTAGAAC",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )


@skipUnless(samtoolsInstalled(), "samtools is not installed")
class TestIgnoreQuality(TestCase, _Mixin):
    """
    Test making majority consensuses with quality ignored.
    """

    ignoreQuality = True


@skipUnless(samtoolsInstalled(), "samtools is not installed")
class TestWithQuality(TestCase, _Mixin):
    """
    Test making majority consensuses with quality considered. The expected
    return reslults for these tests differ from the results that would be
    obtained if quality was ignored, so for that reason these additional tests
    are not part of the common C{_Mixin} class.
    """

    ignoreQuality = False

    def testHighQualityDominates(self):
        """
        If one read has a very high quality base (here 'C', quality ']' = 60)
        that base should take precedence over two other reads that agree with
        each other but which have lower ('5' = 30) quality.
        """
        template = (
            "ACGT",
            "  A",
            "  5",
            "  A",
            "  5",
            "  C",
            "  ]",
        )

        with makeBAM(template) as (fastaFilename, bamFilename):
            self.assertEqual(
                "ACCT",
                consensusFromBAM(
                    bamFilename,
                    quiet=True,
                    referenceFasta=fastaFilename,
                    threshold=0.5,
                    noCoverage="reference",
                    ignoreQuality=self.ignoreQuality,
                ).sequence,
            )
