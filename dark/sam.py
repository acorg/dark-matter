import sys
import numpy as np
from typing import Optional

from json import dump, load
from contextlib import contextmanager
from collections import Counter, defaultdict
from subprocess import CalledProcessError

from dark.utils import pct

from pysam import (
    AlignmentFile,
    AlignedRead,
    CMATCH,
    CINS,
    CDEL,
    CREF_SKIP,
    CSOFT_CLIP,
    CHARD_CLIP,
    CPAD,
    CEQUAL,
    CDIFF,
)

from dark.process import Executor
from dark.fasta import FastaReads
from dark.reads import Read, DNARead


class SamError(Exception):
    "A SAM error."


class UnequalReferenceLengthError(SamError):
    "The references of interest in a SAM/BAM file are not of the same length."


class UnknownReference(SamError):
    "Reference sequence not found in SAM/BAM file."


class UnspecifiedReference(SamError):
    "Reference sequence not specified."


class ReferenceNameMismatchError(SamError):
    "Reference name mismatch."


class InvalidSAM(SamError):
    "SAM/BAM file has unexpected/invalid content."


# From https://samtools.github.io/hts-specs/SAMv1.pdf
CONSUMES_QUERY = {CMATCH, CINS, CSOFT_CLIP, CEQUAL, CDIFF}
CONSUMES_REFERENCE = {CMATCH, CDEL, CREF_SKIP, CEQUAL, CDIFF}


@contextmanager
def samfile(filename):
    """
    A context manager to open and close a SAM/BAM file.

    @param filename: A C{str} file name to open.
    """
    f = AlignmentFile(filename)
    yield f
    f.close()


def samtoolsInstalled():
    """
    Test if samtools is installed.

    @return: A C{bool}, which is C{True} if DIAMOND seems to be installed.
    """
    try:
        Executor().execute("samtools help")
    except CalledProcessError:
        return False
    else:
        return True


def samReferences(filenameOrSamfile):
    """
    List SAM/BAM file reference names.

    @param filenameOrSamfile: Either a C{str} SAM/BAM file name or an
        instance of C{pysam.AlignmentFile}.
    @return: A C{list} of C{str} reference names from the SAM file.
    """

    def _references(sam):
        return [sam.get_reference_name(i) for i in range(sam.nreferences)]

    if isinstance(filenameOrSamfile, str):
        with samfile(filenameOrSamfile) as sam:
            return _references(sam)
    else:
        return _references(filenameOrSamfile)


def samReferencesToStr(filenameOrSamfile, indent=0):
    """
    List SAM/BAM file reference names and lengths.

    @param filenameOrSamfile: Either a C{str} SAM/BAM file name or an
        instance of C{pysam.AlignmentFile}.
    @param indent: An C{int} number of spaces to indent each line.
    @return: A C{str} describing known reference names and their lengths.
    """
    indent = " " * indent

    def _references(sam):
        result = []
        for i in range(sam.nreferences):
            result.append(
                "%s%s (length %d)" % (indent, sam.get_reference_name(i), sam.lengths[i])
            )
        return "\n".join(result)

    if isinstance(filenameOrSamfile, str):
        with samfile(filenameOrSamfile) as sam:
            return _references(sam)
    else:
        return _references(sam)


def _hardClip(sequence, quality, cigartuples):
    """
    Hard clip (if necessary) a sequence.

    @param sequence: A C{str} nucleotide sequence.
    @param quality: A C{str} quality string, or a C{list} of C{int} quality
        values as returned by pysam, or C{None} if the SAM file had a '*'
        for the quality string (which pysam converts to C{None}).
    @param cigartuples: An iterable of (operation, length) tuples, detailing
        the alignment, as per the SAM specification.
    @return: A 3-tuple consisting of
            1) a hard-clipped C{str} sequence if hard-clipping is indicated by
               the CIGAR operations.
            2) a hard-clipped quality C{str} or C{list} (depending on what
               type we were passed) if hard-clipping is indicated by the CIGAR
               operations.
            3) a Boolean, C{True} if hard clipping was performed by this
               function or C{False} if the hard clipping had already been
               done.
    """
    hardClipCount = cigarLength = 0
    for operation, length in cigartuples:
        hardClipCount += operation == CHARD_CLIP
        cigarLength += length if operation in CONSUMES_QUERY else 0

    sequenceLength = len(sequence)
    if quality is not None:
        assert sequenceLength == len(quality)
    clipLeft = clipRight = 0
    clippedSequence = sequence
    clippedQuality = quality

    if sequenceLength > cigarLength:
        alreadyClipped = False
    else:
        assert sequenceLength == cigarLength
        alreadyClipped = True

    if hardClipCount == 0:
        pass
    elif hardClipCount == 1:
        # Hard clip either at the start or the end.
        if cigartuples[0][0] == CHARD_CLIP:
            if not alreadyClipped:
                clipLeft = cigartuples[0][1]
                clippedSequence = sequence[clipLeft:]
                if quality is not None:
                    clippedQuality = quality[clipLeft:]
        elif cigartuples[-1][0] == CHARD_CLIP:
            if not alreadyClipped:
                clipRight = cigartuples[-1][1]
                clippedSequence = sequence[:-clipRight]
                if quality is not None:
                    clippedQuality = quality[:-clipRight]
        else:
            raise ValueError(
                "Invalid CIGAR tuples (%s) contains hard-clipping operation "
                "that is neither at the start nor the end of the sequence."
                % (cigartuples,)
            )
    elif hardClipCount == 2:
        # Hard clip at both the start and end.
        assert cigartuples[0][0] == cigartuples[-1][0] == CHARD_CLIP
        if not alreadyClipped:
            clipLeft, clipRight = cigartuples[0][1], cigartuples[-1][1]
            clippedSequence = sequence[clipLeft:-clipRight]
            if quality is not None:
                clippedQuality = quality[clipLeft:-clipRight]
    else:
        raise ValueError(
            "Invalid CIGAR tuples (%s) specifies hard-clipping %d times (2 "
            "is the maximum)." % (cigartuples, hardClipCount)
        )

    weClipped = bool(clipLeft or clipRight)

    if weClipped:
        assert not alreadyClipped
        if len(clippedSequence) + clipLeft + clipRight != sequenceLength:
            raise ValueError(
                "Sequence %r (length %d) clipped to %r (length %d), but the "
                "difference between these two lengths (%d) is not equal to "
                "the sum (%d) of the left and right clip lengths (%d and %d "
                "respectively). CIGAR tuples: %s"
                % (
                    sequence,
                    len(sequence),
                    clippedSequence,
                    len(clippedSequence),
                    abs(len(sequence) - len(clippedSequence)),
                    clipLeft + clipRight,
                    clipLeft,
                    clipRight,
                    cigartuples,
                )
            )
    else:
        assert len(clippedSequence) == sequenceLength
        if quality is not None:
            assert len(clippedQuality) == sequenceLength

    return clippedSequence, clippedQuality, weClipped


class SAMFilter:
    """
    Filter a SAM/BAM file.

    @param filename: The C{str} name of a BAM/SAM file to filter.
    @param filterRead: A function that takes a C{dark.reads.Read) instance
        and returns either C{None} or a C{Read} instance, according to
        whether the passed read should be omitted or not. If C{None} is
        passed, no filtering on reads (i.e., queries) is done.
    @param referenceIds: Either C{None} or a set of C{str} reference ids
        that should be kept (other references will be dropped).
    @param storeQueryIds: If C{True}, query ids will be stored (in
        C{self.queryIds}) as the SAM/BAM file is read.
    @param dropUnmapped: If C{True}, unmapped matches will be excluded.
    @param dropSecondary: If C{True}, secondary matches will be excluded.
    @param dropSupplementary: If C{True}, supplementary matches will be
        excluded.
    @param dropDuplicates: If C{True}, matches flagged as optical or PCR
        duplicates will be excluded.
    @param keepQCFailures: If C{True}, reads that are considered quality
        control failures will be included.
    @param minScore: If not C{None}, alignments with score tag values less
        than this value will not be output. If given, alignments that do not
        have a score will not be output.
    @param maxScore: If not C{None}, alignments with score tag values greater
        than this value will not be output. If given, alignments that do not
        have a score will not be output.
    @param scoreTag: The alignment tag to extract for minScore and maxScore
        comparisons.
    """

    def __init__(
        self,
        filename,
        filterRead=None,
        referenceIds=None,
        storeQueryIds=False,
        dropUnmapped=False,
        dropSecondary=False,
        dropSupplementary=False,
        dropDuplicates=False,
        keepQCFailures=False,
        minScore=None,
        maxScore=None,
        scoreTag="AS",
    ):
        self.filename = filename
        self.filterRead = filterRead
        self.referenceIds = referenceIds
        self.storeQueryIds = storeQueryIds
        self.dropUnmapped = dropUnmapped
        self.dropSecondary = dropSecondary
        self.dropSupplementary = dropSupplementary
        self.dropDuplicates = dropDuplicates
        self.keepQCFailures = keepQCFailures
        self.minScore = minScore
        self.maxScore = maxScore
        self.scoreTag = scoreTag

        # Detect when there are no filtering criteria, in which case
        # self.filterAlignment can return immediately.
        self.noFiltering = all(
            (
                filterRead is None,
                referenceIds is None,
                not any(
                    (
                        storeQueryIds,
                        dropUnmapped,
                        dropSecondary,
                        dropSupplementary,
                        dropDuplicates,
                        keepQCFailures,
                    )
                ),
                minScore is None,
                maxScore is None,
            )
        )

    @staticmethod
    def addFilteringOptions(
        parser,
        samfileIsPositional=False,
        samfileAction="store",
        samfileRequired=True,
        samfileNargs=None,
        referenceIdRequired=False,
    ):
        """
        Add options to an argument parser for filtering SAM/BAM.

        @param parser: An C{argparse.ArgumentParser} instance.
        @param samfileIsPositional: If C{True} the SAM/BAM file must
            be given as the final argument on the command line (without
            being preceded by --samfile).
        @param samfileAction: A C{str} action to take when 'samfile' arguments
            are found on the command line. Pass 'append' to allow multiple SAM
            files.
        @param samfileRequired: If C{True}, --samfile must be given on the
            command line. This is only relevant when samfileIsPositional is
            C{False} because positional arguments are always required. It may
            seem strange to allow no --samfile argument in a class that filters
            SAM files, but this option can be used in conjuction with scripts
            that can optionally take a SAM file (e.g.,
            genome-protein-summary.py).
        @param samfileNargs: The value to pass for 'nargs' in adding the
            samfile option.
        @param referenceIdRequired: If C{True}, make the --referenceId option
            required.
        """
        if samfileIsPositional:
            # Positional arguments are always required.
            assert (
                samfileRequired
            ), "samfileIsPositional is True, so samfileRequired must also be True."
            if samfileNargs is None:
                parser.add_argument(
                    "samfile", action=samfileAction, help="The SAM/BAM file to filter."
                )
            else:
                parser.add_argument(
                    "samfile",
                    action=samfileAction,
                    nargs=samfileNargs,
                    help="The SAM/BAM file to filter.",
                )
        else:
            if samfileNargs is None:
                parser.add_argument(
                    "--samfile",
                    required=samfileRequired,
                    action=samfileAction,
                    help="The SAM/BAM file to filter.",
                )
            else:
                parser.add_argument(
                    "--samfile",
                    required=samfileRequired,
                    nargs=samfileNargs,
                    action=samfileAction,
                    help="The SAM/BAM file to filter.",
                )

        parser.add_argument(
            "--referenceId",
            metavar="ID",
            action="append",
            required=referenceIdRequired,
            help=(
                "A reference sequence id whose alignments should be kept "
                "(alignments against other references will be dropped). "
                "If omitted, alignments against all references will be "
                "kept. May be repeated."
            ),
        )

        parser.add_argument(
            "--dropUnmapped",
            action="store_true",
            help="If given, unmapped matches will not be output.",
        )

        parser.add_argument(
            "--dropSecondary",
            action="store_true",
            help="If given, secondary matches will not be output.",
        )

        parser.add_argument(
            "--dropSupplementary",
            action="store_true",
            help="If given, supplementary matches will not be output.",
        )

        parser.add_argument(
            "--dropDuplicates",
            action="store_true",
            help=(
                "If given, matches flagged as optical or PCR duplicates "
                "will not be output."
            ),
        )

        parser.add_argument(
            "--keepQCFailures",
            action="store_true",
            help=(
                "If given, reads that are considered quality control "
                "failures will be included in the output."
            ),
        )

        parser.add_argument(
            "--minScore",
            type=float,
            metavar="FLOAT",
            help=(
                "If given, alignments with --scoreTag (default AS) values "
                "less than this value will not be output. If given, "
                "alignments that do not have a score will not be output."
            ),
        )

        parser.add_argument(
            "--maxScore",
            type=float,
            metavar="FLOAT",
            help=(
                "If given, alignments with --scoreTag (default AS) values "
                "greater than this value will not be output. If given, "
                "alignments that do not have a score will not be output."
            ),
        )

        parser.add_argument(
            "--scoreTag",
            default="AS",
            metavar="TAG",
            help=(
                "The alignment tag to extract for --minScore and --maxScore "
                "comparisons."
            ),
        )

    @classmethod
    def parseFilteringOptions(cls, args, filterRead=None, storeQueryIds=False):
        """
        Parse command line options (added in C{addSAMFilteringOptions}.

        @param args: The command line arguments, as returned by
            C{argparse.parse_args}.
        @param filterRead: A one-argument function that accepts a read
            and returns C{None} if the read should be omitted in filtering
            or else a C{Read} instance.
        @param storeQueryIds: If C{True}, query ids will be stored as the
            SAM/BAM file is read.
        @return: A C{SAMFilter} instance.
        """
        return cls(
            args.samfile,
            filterRead=filterRead,
            referenceIds=set(args.referenceId) if args.referenceId else None,
            storeQueryIds=storeQueryIds,
            dropUnmapped=args.dropUnmapped,
            dropSecondary=args.dropSecondary,
            dropSupplementary=args.dropSupplementary,
            dropDuplicates=args.dropDuplicates,
            keepQCFailures=args.keepQCFailures,
            minScore=args.minScore,
            maxScore=args.maxScore,
        )

    def filterAlignment(self, alignment):
        """
        Test an alignment to see if it passes all filtering criteria.

        @param alignment: A pysam alignment instance.
        @return: A C{bool}, C{True} if the alignment passes our filtering,
            C{False} if it should be discarded.
        """
        if self.noFiltering:
            return True

        if self.minScore is not None or self.maxScore is not None:
            try:
                score = alignment.get_tag(self.scoreTag)
            except KeyError:
                return False
            else:
                if (self.minScore is not None and score < self.minScore) or (
                    self.maxScore is not None and score > self.maxScore
                ):
                    return False

        return (
            self.filterRead is None
            or self.filterRead(
                Read(alignment.query_name, alignment.query_sequence, alignment.qual)
            )
        ) and not (
            (self.referenceIds and alignment.reference_name not in self.referenceIds)
            or (alignment.is_unmapped and self.dropUnmapped)
            or (alignment.is_secondary and self.dropSecondary)
            or (alignment.is_supplementary and self.dropSupplementary)
            or (alignment.is_duplicate and self.dropDuplicates)
            or (alignment.is_qcfail and not self.keepQCFailures)
        )

    def alignments(self):
        """
        Get alignments from the SAM/BAM file, subject to filtering.

        @return: A generator that yields pysam alignment instances that pass
            our filtering criteria.
        """
        if self.storeQueryIds:
            queryIds = self.queryIds = set()
        else:
            queryIds = None

        lastAlignment = None
        count = 0
        with samfile(self.filename) as samAlignment:
            for count, alignment in enumerate(samAlignment.fetch(), start=1):
                if queryIds is not None:
                    queryIds.add(alignment.query_name)

                # Secondary and supplementary alignments may have a '*'
                # (pysam returns this as None) SEQ field, indicating that
                # the previous sequence should be used. This is best
                # practice according to section 2.5.2 of
                # https://samtools.github.io/hts-specs/SAMv1.pdf So we use
                # the last alignment query and quality strings if we get
                # None as a query sequence.
                if alignment.query_sequence is None:
                    if lastAlignment is None:
                        raise InvalidSAM(
                            "pysam produced an alignment (number %d) with no "
                            "query sequence without previously giving an "
                            "alignment with a sequence." % count
                        )
                    # Use the previous query sequence and quality. I'm not
                    # making the call to _hardClip dependent on
                    # alignment.cigartuples (as in the else clause below)
                    # because I don't think it's possible for
                    # alignment.cigartuples to be None in this case. If we
                    # have a second match on a query, then it must be
                    # aligned to something (i.e., it cannot be unmapped
                    # with no CIGAR string). The assertion will tell us if
                    # this is ever not the case.
                    assert alignment.cigartuples
                    (
                        alignment.query_sequence,
                        alignment.query_qualities,
                        _,
                    ) = _hardClip(
                        lastAlignment.query_sequence,
                        lastAlignment.query_qualities,
                        alignment.cigartuples,
                    )
                else:
                    lastAlignment = alignment
                    if alignment.cigartuples:
                        (
                            alignment.query_sequence,
                            alignment.query_qualities,
                            _,
                        ) = _hardClip(
                            alignment.query_sequence,
                            alignment.query_qualities,
                            alignment.cigartuples,
                        )

                if self.filterAlignment(alignment):
                    yield alignment

        self.alignmentCount = count

    def referenceLengths(self):
        """
        Get the lengths of wanted references.

        @raise UnknownReference: If a reference id is not present in the
            SAM/BAM file.
        @return: A C{dict} of C{str} reference id to C{int} length with a key
            for each reference id in C{self.referenceIds} or for all references
            if C{self.referenceIds} is C{None}.
        """
        result = {}
        with samfile(self.filename) as sam:
            if self.referenceIds:
                for referenceId in self.referenceIds:
                    tid = sam.get_tid(referenceId)
                    if tid == -1:
                        raise UnknownReference(
                            "Reference %r is not present in the SAM/BAM file."
                            % referenceId
                        )
                    else:
                        result[referenceId] = sam.lengths[tid]
            else:
                result = dict(zip(sam.references, sam.lengths))

        return result


class PaddedSAM:
    """
    Obtain aligned (padded) queries from a SAM/BAM file.

    @param samFilter: A C{SAMFilter} instance.
    @raises UnequalReferenceLengthError: If C{referenceName} is C{None}
        and the reference sequence lengths in the SAM/BAM file are not all
        identical.
    @raises UnknownReference: If C{referenceName} does not exist.
    """

    def __init__(self, samFilter):
        referenceLengths = samFilter.referenceLengths()

        if len(set(referenceLengths.values())) != 1:
            raise UnequalReferenceLengthError(
                "Your %d SAM/BAM file reference sequence "
                "lengths (%s) are not all identical."
                % (
                    len(referenceLengths),
                    ", ".join(
                        "%s=%d" % (id_, referenceLengths[id_])
                        for id_ in sorted(referenceLengths)
                    ),
                )
            )

        # Get the length of any of the sequences (they are all identical).
        self.referenceLength = referenceLengths.popitem()[1]
        self.samFilter = samFilter
        # self.referenceInsertions will be keyed by query id (the query
        # that would cause a reference insertion). The values will be lists
        # of 2-tuples, with each 2-tuple containing an offset into the
        # reference sequence and the C{str} of nucleotides that would be
        # inserted starting at that offset.
        self.referenceInsertions = defaultdict(list)

    def queries(
        self,
        rcSuffix="",
        rcNeeded=False,
        padChar="-",
        queryInsertionChar="N",
        unknownQualityChar="!",
        allowDuplicateIds=False,
        addAlignment=False,
    ):
        """
        Produce padded (with gaps) queries according to the CIGAR string and
        reference sequence length for each matching query sequence.

        @param rcSuffix: A C{str} to add to the end of query names that are
            reverse complemented. This is added before the /1, /2, etc., that
            are added for duplicated ids (if there are duplicates and
            C{allowDuplicateIds} is C{False}.
        @param rcNeeded: If C{True}, queries that are flagged as matching when
            reverse complemented should have reverse complementing when
            preparing the output sequences. This must be used if the program
            that created the SAM/BAM input flags reversed matches but does not
            also store the reverse complemented query.
        @param padChar: A C{str} of length one to use to pad queries with to
            make them the same length as the reference sequence.
        @param queryInsertionChar:  A C{str} of length one to use to insert
            into queries when the CIGAR string indicates that the alignment
            of a query would cause a deletion in the reference. This character
            is inserted as a 'missing' query character (i.e., a base that can
            be assumed to have been lost due to an error) whose existence is
            necessary for the match to continue.
        @param unknownQualityChar: The character to put into the quality
            string when unknown bases are inserted in the query or the query
            is padded on the left/right with gaps.
        @param allowDuplicateIds: If C{True}, repeated query ids (due to
            secondary or supplemental matches) will not have /1, /2, etc.
            appended to their ids. So repeated ids may appear in the yielded
            FASTA.
        @param addAlignment: If C{True} the reads yielded by the returned
            generator will also have an C{alignment} attribute, being the
            C{pysam.AlignedSegment} for the query.
        @raises InvalidSAM: If a query has an empty SEQ field and either there
            is no previous alignment or the alignment is not marked as
            secondary or supplementary.
        @return: A generator that yields C{Read} instances that are padded
            with gap characters to align them to the length of the reference
            sequence. See C{addAlignment}, above, to yield reads with the
            corresponding C{pysam.AlignedSegment}.
        """
        referenceLength = self.referenceLength

        # Hold the count for each id so we can add /1, /2 etc to duplicate
        # ids (unless --allowDuplicateIds was given).
        idCount = Counter()

        MATCH_OPERATIONS = {CMATCH, CEQUAL, CDIFF}

        for lineNumber, alignment in enumerate(self.samFilter.alignments(), start=1):
            query = alignment.query_sequence
            quality = "".join(chr(q + 33) for q in alignment.query_qualities)

            if alignment.is_reverse:
                if rcNeeded:
                    query = DNARead("id", query).reverseComplement().sequence
                    quality = quality[::-1]
                if rcSuffix:
                    alignment.query_name += rcSuffix

            # Adjust the query id if it's a duplicate and we're not allowing
            # duplicates.
            if allowDuplicateIds:
                queryId = alignment.query_name
            else:
                count = idCount[alignment.query_name]
                idCount[alignment.query_name] += 1
                queryId = alignment.query_name + ("" if count == 0 else "/%d" % count)

            referenceStart = alignment.reference_start
            atStart = True
            queryIndex = 0
            referenceIndex = referenceStart
            alignedSequence = ""
            alignedQuality = ""

            for operation, length in alignment.cigartuples:
                # The operations are tested in the order they appear in
                # https://samtools.github.io/hts-specs/SAMv1.pdf It would be
                # more efficient to test them in order of frequency of
                # occurrence.
                if operation in MATCH_OPERATIONS:
                    atStart = False
                    alignedSequence += query[queryIndex : queryIndex + length]
                    alignedQuality += quality[queryIndex : queryIndex + length]
                elif operation == CINS:
                    # Insertion to the reference. This consumes query bases but
                    # we don't output them because the reference cannot be
                    # changed.  I.e., these bases in the query would need to be
                    # inserted into the reference.  Remove these bases from the
                    # query but record what would have been inserted into the
                    # reference.
                    atStart = False
                    self.referenceInsertions[queryId].append(
                        (referenceIndex, query[queryIndex : queryIndex + length])
                    )
                elif operation == CDEL:
                    # Delete from the reference. Some bases from the reference
                    # would need to be deleted to continue the match. So we put
                    # an insertion into the query to compensate.
                    atStart = False
                    alignedSequence += queryInsertionChar * length
                    alignedQuality += unknownQualityChar * length
                elif operation == CREF_SKIP:
                    # Skipped reference. Opens a gap in the query. For
                    # mRNA-to-genome alignment, an N operation represents an
                    # intron.  For other types of alignments, the
                    # interpretation of N is not defined. So this is unlikely
                    # to occur.
                    atStart = False
                    alignedSequence += queryInsertionChar * length
                    alignedQuality += unknownQualityChar * length
                elif operation == CSOFT_CLIP:
                    # Bases in the query that are not part of the match. We
                    # remove these from the query if they protrude before the
                    # start or after the end of the reference. According to the
                    # SAM docs, 'S' operations may only have 'H' operations
                    # between them and the ends of the CIGAR string.
                    if atStart:
                        # Don't set atStart=False, in case there's another 'S'
                        # operation.
                        unwantedLeft = length - referenceStart
                        if unwantedLeft > 0:
                            # The query protrudes left. Copy its right part.
                            alignedSequence += query[
                                queryIndex + unwantedLeft : queryIndex + length
                            ]
                            alignedQuality += quality[
                                queryIndex + unwantedLeft : queryIndex + length
                            ]
                            referenceStart = 0
                        else:
                            referenceStart -= length
                            alignedSequence += query[queryIndex : queryIndex + length]
                            alignedQuality += quality[queryIndex : queryIndex + length]
                    else:
                        unwantedRight = (
                            referenceStart + len(alignedSequence) + length
                        ) - referenceLength

                        if unwantedRight > 0:
                            # The query protrudes right. Copy its left part.
                            alignedSequence += query[
                                queryIndex : queryIndex + length - unwantedRight
                            ]
                            alignedQuality += quality[
                                queryIndex : queryIndex + length - unwantedRight
                            ]
                        else:
                            alignedSequence += query[queryIndex : queryIndex + length]
                            alignedQuality += quality[queryIndex : queryIndex + length]
                elif operation == CHARD_CLIP:
                    # Some bases have been completely removed from the query.
                    # This (H) can only be present as the first and/or last
                    # operation. There is nothing to do as the bases are simply
                    # not present in the query string in the SAM/BAM file.
                    pass
                elif operation == CPAD:
                    # This is "silent deletion from the padded reference",
                    # which consumes neither query nor reference.
                    atStart = False
                else:
                    raise ValueError("Unknown CIGAR operation:", operation)

                if operation in CONSUMES_QUERY:
                    queryIndex += length

                if operation in CONSUMES_REFERENCE:
                    referenceIndex += length

            if queryIndex != len(query):
                # Oops, we did not consume the entire query.
                raise ValueError(
                    "Query %r not fully consumed when parsing CIGAR string. "
                    "Query %r (len %d), final query index %d, CIGAR: %r"
                    % (
                        alignment.query_name,
                        query,
                        len(query),
                        queryIndex,
                        alignment.cigartuples,
                    )
                )

            # We cannot test we consumed the entire reference.  The CIGAR
            # string applies to (and exhausts) the query but is silent
            # about the part of the reference that lies to the right of the
            # aligned query.

            # Put gap characters before and after the aligned sequence so that
            # it is offset properly and matches the length of the reference.
            padRightLength = referenceLength - (referenceStart + len(alignedSequence))
            paddedSequence = (
                padChar * referenceStart + alignedSequence + padChar * padRightLength
            )
            paddedQuality = (
                unknownQualityChar * referenceStart
                + alignedQuality
                + unknownQualityChar * padRightLength
            )

            read = Read(queryId, paddedSequence, paddedQuality)

            if addAlignment:
                read.alignment = alignment

            yield read


class DistanceMatrix:
    """
    Maintain score (matching) information between reads and reference sequences
    as read from a SAM file (or files) and provide methods for calculating
    distance between references.
    """

    def __init__(self):
        # self.scores is a similarity matrix between references and
        # queries from a SAM file. The values are 1.0 or 0.0 when queries
        # match / do not match references unless we are using alignment
        # scores, as generated by tools like bowtie2 and bwa and stored in
        # a SAM tag. In that case, the alignment score value is stored (see
        # https://en.wikipedia.org/wiki/SAM_(file_format)#Optional_fields
        # for SAM tag info).
        #
        # The dictionary is keyed by str reference ids and the values are
        # dicts keyed by str query ids. The values of the sub-dictionaries
        # are the float scores (similarities). E.g.
        #
        # {
        #     'reference-1': {
        #         'query-1': 77.5,
        #         'query-2': 37.0,
        #     },
        #     'reference-2': {
        #         'query-1': 78.0,
        #     },
        # }
        #
        # If a reference/query id pair is not present in self.scores, it
        # means that query did not match that reference. The score method
        # (below) will assign such a pair a similarity score of 0.0.
        self.scores = {}

    def addFile(self, filename, scoreTag=None):
        """
        Add information from a SAM file.

        @param filename: A C{str} SAM file name to read.
        @param scoreTag: A C{str} alignment score tag to use instead of binary
            match / non-match values that are used if C{scoreTag} is C{None}.
            The score tag must be present in the optional SAM fields.
        """
        if scoreTag:

            def getScore(alignment, filename, count):
                """
                Get the alignment score from a SAM alignment.

                @param alignment: A samtools alignment.
                @param filename: The C{str} SAM file name being read.
                @param count: The C{int} line number in the SAM file.
                """
                try:
                    score = alignment.get_tag(scoreTag)
                except KeyError:
                    raise ValueError(
                        f"Alignment {count} in {filename!r} has no "
                        f"{scoreTag!r} score tag."
                    )
                else:
                    if score < 0.0:
                        raise ValueError(
                            f"Alignment {count} in {filename!r} has tag "
                            f"{scoreTag!r} with negative value ({score})."
                        )
                    return score

        else:

            def getScore(alignment, filename, count):
                """
                Return a binary score of 1.0 seeing as we know the query was
                mapped and no score tag is in use.

                @param args: ignored.
                """
                return 1.0

        scores = self.scores

        with samfile(filename) as samAlignment:
            for count, alignment in enumerate(samAlignment.fetch(), start=1):
                if alignment.is_unmapped:
                    continue

                queryId = alignment.query_name
                referenceId = alignment.reference_name
                score = getScore(alignment, filename, count)

                try:
                    preExisting = scores[referenceId][queryId]
                except KeyError:
                    if referenceId not in scores:
                        scores[referenceId] = {}

                    scores[referenceId][queryId] = score
                else:
                    if score > preExisting:
                        scores[referenceId][queryId] = score

    def score(self, referenceId, queryId):
        """
        Get the score for a reference and query.

        Although this method name does not start with an underscore, it's not
        likely to be used directly when this class is used to calculate
        distances between references according to query matching. The scores
        between references and queries (returned by this function) are just
        a component of the overall distance.

        @param referenceId: A C{str} reference id.
        @param queryId: A C{str} query id.
        @return: A C{float} score (zero if there was no match between the query
            and reference) or if the reference id or query id have never been
            seen. This could arguably cause a KeyError, but I decided against
            that to keep the code simpler / faster. If you pass in unknown
            query or reference ids, that's your problem (this should be very
            unlikely because our caller is operating on a SAM file, so the
            place they're most likely to be getting query and reference ids
            from to pass us is from that same SAM file.
        """
        try:
            return self.scores[referenceId][queryId]
        except KeyError:
            return 0.0

    def jaccardDistance(self, referenceId1, referenceId2):
        """
        Get the Jaccard distance between two references.  See
        https://en.wikipedia.org/wiki/Jaccard_index

        @param referenceId1: A C{str} reference id.
        @param referenceId2: A C{str} reference id.
        @return: A C{float} distance.
        """
        queryIds1 = set(self.scores.get(referenceId1, []))
        queryIds2 = set(self.scores.get(referenceId2, []))

        denominator = len(queryIds1 | queryIds2)

        if denominator:
            numerator = len(queryIds1 & queryIds2)
            similarity = numerator / denominator
            distance = 1.0 - similarity
            return distance
        else:
            # The references have no matching queries in common.
            return 1.0

    def soergelDistance(self, referenceId1, referenceId2):
        """
        Get the Soergel (i.e., weighted Jaccard) distance between two
        references. See https://en.wikipedia.org/wiki/Jaccard_index

        Note that this will be the same as the jaccardDistance unless a
        scoreTag was used in reading the SAM file.

        @param referenceId1: A C{str} reference id.
        @param referenceId2: A C{str} reference id.
        @return: A C{float} distance.
        """
        queryIds1 = set(self.scores.get(referenceId1, []))
        queryIds2 = set(self.scores.get(referenceId2, []))

        score = self.score
        scores = [
            (score(referenceId1, queryId), score(referenceId2, queryId))
            for queryId in (queryIds1 | queryIds2)
        ]

        denominator = sum(max(a, b) for a, b in scores)

        if denominator:
            numerator = sum(min(a, b) for a, b in scores)
            similarity = numerator / denominator
            distance = 1.0 - similarity
            return distance
        else:
            # Neither reference had a positive match against any query.
            return 1.0

    def matrix(
        self, referenceIds=None, metric="soergel", similarity=False, returnDict=False
    ):
        """
        Compute a distance matrix.

        @param referenceIds: An iterable of C{str} reference ids. If C{None},
            all known reference ids will be used.
        @param metric: A C{str}, either 'soergel' or 'jaccard'.
        @param similarity: If C{True}, return a similarity matrix.
        @param returnDict: If C{True}, return a C{dict} (see @return, below).
        @return: If C{returnDict} is C{False}, return a symmetric square
            C{np.array} of C{float} inter-reference distances (or similarities
            if C{similarity} is C{True}) in the range (0.0, 1.0). The order of
            the entries in the rows and columns of the returned array matches
            that of C{referenceIds} (if given) or the order in C{self.scores}
            otherwise.  If C{returnDict} is C{True}, return a C{dict} of
            C{dict}s, indexed by the two reference ids, with values as above.
        """
        assert metric in {"jaccard", "soergel"}
        referenceIds = tuple(referenceIds or self.scores)
        nIds = len(referenceIds)
        matrix = defaultdict(dict) if returnDict else np.empty((nIds, nIds))
        identityValue = 1.0 if similarity else 0.0
        func = self.jaccardDistance if metric == "jaccard" else self.soergelDistance

        if similarity:

            def metricFunc(ref1, ref2):
                return 1.0 - func(ref1, ref2)

        else:

            def metricFunc(ref1, ref2):
                return func(ref1, ref2)

        for index1, ref1 in enumerate(referenceIds):
            for index2, ref2 in enumerate(referenceIds):
                i, j = (ref1, ref2) if returnDict else (index1, index2)
                value = identityValue if i == j else metricFunc(ref1, ref2)
                matrix[i][j] = matrix[j][i] = value

        return matrix

    def save(self, fp):
        """
        Save the similarity matrix as JSON.

        @param fp: An open file pointer to write to.
        """
        dump(self.scores, fp, sort_keys=True, indent=4)

    def load(self, fp):
        """
        Load a similarity matrix from a JSON file.

        @param fp: An open file pointer to read from.
        """
        self.scores = load(fp)


def getReferenceInfo(
    bam, bamFilename, bamId=None, referenceFasta=None, fastaId=None, quiet=False
):
    """
    Get the id and length of a BAM file reference sequence.

    @param bam: An open BAM file.
    @param bamFilename: The C{str} file name of the BAM file.
    @param bamId: A C{str} BAM file reference name indicating which aligned
        reads to make a consensus from. If not given, will be inferred
        from the BAM file header.
    @param referenceFasta: A C{str} file name containing the sequence that was
        aligned to in making the BAM file.
    @param fastaId: A C{str} reference name indicating which sequence in
        C{referenceFasta} to use as a reference. Only considered if
        C{referenceFasta} is given. If not given and C{referenceFasta} is,
        the reference id will be inferred from reference names in the BAM
        header, or will be taken as the id of the first sequence in
        C{referenceFasta}.
    @param quiet: If C{True}, suppress diagnostic output.
    @raise UnknownReference: For an unknown reference.
    @raise UnspecifiedReference: If a reference is not given and cannot be
        inferred.
    @raise ReferenceNameMismatchError: If the names of the reference in the BAM
        file and the FASTA reference file (if given) do not match (unless both
        are explicitly given, in which case mismatching names are allowed as
        long as the sequence lengths are the same).
    @raise UnequalReferenceLengthError: If the lengths of the reference in the
        BAM file and the FASTA reference file (if given) do not match
    @return: A 2-tuple of C{str}, C{int} giving the refrence sequence id and
        length.
    """
    bamReferences = set(samReferences(bam))
    checkName = checkLength = True
    inferredBamId = bamId is None
    inferredFastaId = fastaId is None
    referenceLength = None

    if bamId is None:
        inferredBamId = True
    else:
        if bamId not in bamReferences:
            raise UnknownReference(
                f"BAM file {str(bamFilename)!r} does not mention a reference "
                f"with id {bamId!r}. Known references are: "
                f'{", ".join(sorted(bamReferences))}.'
            )

        inferredBamId = False
        tid = bam.get_tid(bamId)
        referenceLength = bam.lengths[tid]

        if fastaId is not None:
            # Both ids have been given explicitly, so don't check that they
            # are identical, assume the user knows what they are doing (we
            # will just check the lengths, if necessary).
            checkName = False

    if referenceFasta is None:
        reference = None
    else:
        if fastaId is None:
            # We're given a reference FASTA file, but no id to use. Look at
            # the sequences in the FASTA and take the first one that has a
            # name and length matching a BAM reference (we could actually
            # look to see if there are more than one and, if so, give an
            # error if their sequences are not identical).
            firstRead = None
            for read in FastaReads(referenceFasta):
                if firstRead is None:
                    firstRead = read
                if read.id in bamReferences:
                    tid = bam.get_tid(read.id)
                    if len(read) == bam.lengths[tid]:
                        reference = read
                        checkLength = checkName = False
                        break
            else:
                # Try the first sequence in the FASTA file.
                reference = firstRead

            if firstRead is None:
                raise UnknownReference(
                    f"The FASTA reference file {str(referenceFasta)!r} "
                    f"contained no sequences."
                )
        else:
            for read in FastaReads(referenceFasta):
                if read.id == fastaId:
                    reference = read
                    break
            else:
                raise UnknownReference(
                    f"No sequence with id {fastaId!r} "
                    f"found in {str(referenceFasta)!r}."
                )

        # Set reference length if it has not already been set due to a
        # given bamId.
        if referenceLength is None:
            referenceLength = len(reference)

    if bamId is None:
        if len(bamReferences) == 1:
            # This is the only possibility.
            bamId = tuple(bamReferences)[0]
            tid = bam.get_tid(bamId)
            referenceLength = bam.lengths[tid]
        elif reference and reference.id in bamReferences:
            tid = bam.get_tid(reference.id)
            if len(reference) == bam.lengths[tid]:
                bamId = reference.id
                checkLength = checkName = False
        else:
            raise UnspecifiedReference(
                f"Could not infer a BAM reference. Available references are: "
                f'{", ".join(sorted(bamReferences))}.'
            )

    if inferredBamId and not quiet:
        print(f"BAM reference id {bamId!r} inferred from context.", file=sys.stderr)

    if reference:
        if inferredFastaId and not quiet:
            print(
                f"FASTA reference id {reference.id!r} inferred from " f"context.",
                file=sys.stderr,
            )

        if checkName and reference.id != bamId:
            raise ReferenceNameMismatchError(
                f"Reference FASTA sequence name {reference.id!r} does "
                f"not match the BAM id {bamId!r}."
            )

        if checkLength and len(reference) != referenceLength:
            raise UnequalReferenceLengthError(
                f"Reference FASTA sequence {reference.id!r} has length "
                f"{len(reference)}, but the BAM reference {bamId!r} has "
                f"length {referenceLength}."
            )

    return bamId, reference, referenceLength


class ReferenceReads:
    """
    Hold information about the reads that match a reference.

    @param id_: A C{str} reference id.
    @param length: The C{int} length of the reference.
    """

    def __init__(self, id_: str, length: int) -> None:
        assert length
        self.id_ = id_
        self.length = length
        self.coveredOffsets: set = set()
        self.duplicate: set = set()
        self.nonDuplicate: set = set()
        self.primary: set = set()
        self.qcFail: set = set()
        self.readIds: set = set()
        self.secondary: set = set()
        self.supplementary: set = set()

    def add(self, read: AlignedRead) -> None:
        """
        Add a matching read to this reference.

        @param read: An C{AlignedRead} to add.
        """
        id_ = read.query_name

        if read.is_secondary:
            self.secondary.add(id_)
        elif read.is_supplementary:
            self.supplementary.add(id_)
        else:
            self.primary.add(id_)

        # Note that it is possible that a SAM/BAM file has the same read id
        # more than once with identical flags (with the flags not including
        # 0x400 (1024, pysam.FDUP). Bowtie2 produces such lines. In such
        # cases, the read has not been flagged as a duplicate even though
        # it clearly is. To deal with this, if we see a read id again, we
        # unconditionally put it into self.duplicate (and make sure it is
        # not in self.nonDuplicate).
        if id_ in self.readIds:
            self.duplicate.add(id_)
            self.nonDuplicate.discard(id_)
        else:
            if read.is_duplicate:
                self.duplicate.add(id_)
            else:
                self.nonDuplicate.add(id_)

        self.readIds.add(id_)

        if read.is_qcfail:
            self.qcFail.add(id_)

        if read.reference_length is not None:
            self.coveredOffsets.update(
                range(
                    read.reference_start,
                    read.reference_start + read.reference_length,
                )
            )

    def coverage(self) -> float:
        """
        Determine the fraction of the reference covered by at least one read.

        @return: The C{float} coverage fraction.
        """
        return len(self.coveredOffsets) / self.length


class ReadsSummary:
    """
    Hold information about reads found in a BAM/SAM file.

    @param samFile: A C{str} filename to read.
    """

    def __init__(self, samFile: str) -> None:
        self.mapped: set[str] = set()
        self.unmapped: set[str] = set()
        self.readIds: set[str] = set()
        self.references: dict[str, ReferenceReads] = {}
        self.referenceLengths = SAMFilter(samFile).referenceLengths()
        self.sortedBy: Optional[str] = None
        self.sortedReferences: list[ReferenceReads] = []

        with samfile(samFile) as fp:
            for read in fp.fetch():
                id_ = read.query_name
                self.readIds.add(id_)

                if read.is_unmapped:
                    self.unmapped.add(id_)
                    continue

                self.mapped.add(id_)

                referenceName = read.reference_name
                try:
                    reference = self.references[referenceName]
                except KeyError:
                    reference = self.references[referenceName] = ReferenceReads(
                        referenceName, self.referenceLengths[referenceName]
                    )

                reference.add(read)

    def toString(
        self,
        excludeZeroes: bool = False,
        excludeIfNoAdditional: bool = False,
        sortBy: str = "count",
    ) -> str:
        """
        Return a string suitable for printing.

        @param excludeZeroes: Do not summarize references unmatched by any read.
        @param excludeIfNoAdditional: Do not summarize references that are mapped
            by no additional reads (i.e., reads other than those matching
            already-summarized references). This is useful C{sortBy} is 'count',
            making it possible to not summarize references only matched by reads
            matching references that have already been summarized
        @param sortBy: A C{str} indicating how to sort. Use 'count' to sort by
            decreasing read count (i.e., the reference with the highest number
            of matching reads comes first), or 'coverage' to sort by decreasing
            genome coverage, or 'name' to sort by reference name.
        @return: A C{list} of sorted C{ReferenceReads} instances.
        """
        totalReads = len(self.readIds)
        result = [
            "Found a total of %d read%s (%d mapped, %d unmapped) mapping "
            "against %d of %d reference%s."
            % (
                totalReads,
                "" if totalReads == 1 else "s",
                len(self.mapped),
                len(self.unmapped),
                len(self.references),
                len(self.referenceLengths),
                "" if len(self.referenceLengths) == 1 else "s",
            )
        ]

        if not self.references:
            result.append("No reads were mapped to any reference.")
            return "\n".join(result)

        cumulativeReadIds: set[str] = set()

        if self.sortedBy != sortBy:
            self.sortedReferences = self.sortReferences(sortBy)

        for count, reference in enumerate(self.sortedReferences, start=1):
            readIds = reference.readIds
            readCount = len(readIds)
            if readCount == 0 and excludeZeroes:
                continue
            newReadCount = len(readIds - cumulativeReadIds)
            if newReadCount == 0 and excludeIfNoAdditional:
                continue
            cumulativeReadIds.update(readIds)
            result.append(
                "\nReference %d: %s (%d nt):\n"
                "  Genome coverage: %.2f%%\n"
                "  Overall reads mapped to the reference: %s\n"
                "    Non-duplicates: %s, Duplicates: %s, QC fails: %s\n"
                "    Primary: %s, Secondary: %s, Supplementary: %s\n"
                "    Reads not matching any reference above: %s\n"
                "  Previously unmatched reads for this reference: %s"
                % (
                    count,
                    reference.id_,
                    reference.length,
                    reference.coverage() * 100.0,
                    pct(readCount, totalReads),
                    pct(len(reference.nonDuplicate), readCount),
                    pct(len(reference.duplicate), readCount),
                    pct(len(reference.qcFail), readCount),
                    pct(len(reference.primary), readCount),
                    pct(len(reference.secondary), readCount),
                    pct(len(reference.supplementary), readCount),
                    pct(newReadCount, totalReads),
                    pct(newReadCount, readCount),
                )
            )

        return "\n".join(result)

    def sortReferences(self, sortBy: str = "count"):
        """
        Sort the references.

        @param sortBy: A C{str} indicating how to sort. Use 'count' to sort by
            decreasing read count (i.e., the reference with the highest number
            of matching reads comes first), or 'coverage' to sort by decreasing
            genome coverage, or 'name' to sort by reference name.
        @raise AssertionError: If C{sortBy} has an unknown value.
        @return: A C{list} of sorted C{ReferenceReads} instances.
        """
        if self.sortedBy == sortBy:
            return self.sortedReferences

        if sortBy == "count":
            reverse = True

            def key(reference):
                return len(reference.readIds), reference.coverage()

        elif sortBy == "coverage":
            reverse = True

            def key(reference):
                return reference.coverage(), len(reference.readIds)

        elif sortBy == "name":
            reverse = False

            def key(reference):
                return reference.id_, len(reference.readIds), reference.coverage()

        else:
            raise ValueError(f"Unknown sortBy value {sortBy!r}.")

        self.sortedReferences = sorted(
            self.references.values(), key=key, reverse=reverse
        )
        self.sortedBy = sortBy

        return self.sortedReferences
