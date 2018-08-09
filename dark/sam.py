import six

from itertools import chain
from contextlib import contextmanager
from collections import Counter, defaultdict

from pysam import (
    AlignmentFile, CMATCH, CINS, CDEL, CREF_SKIP, CSOFT_CLIP, CHARD_CLIP, CPAD,
    CEQUAL, CDIFF)

from dark.reads import Read, DNARead


class UnequalReferenceLengthError(Exception):
    "The references of interest in a SAM/BAM file are not of the same length."


class UnknownReference(Exception):
    "Reference sequence not found in SAM/BAM file."


class InvalidSAM(Exception):
    "SAM/BAM file has unexpected/invalid content."


# From https://samtools.github.io/hts-specs/SAMv1.pdf
_CONSUMES_QUERY = {CMATCH, CINS, CSOFT_CLIP, CEQUAL, CDIFF}
_CONSUMES_REFERENCE = {CMATCH, CDEL, CREF_SKIP, CEQUAL, CDIFF}


@contextmanager
def samfile(filename):
    """
    A context manager to open and close a SAM/BAM file.

    @param filename: A C{str} file name to open.
    """
    f = AlignmentFile(filename)
    yield f
    f.close()


def samReferencesToStr(filenameOrSamfile, indent=0):
    """
    List SAM/BAM file reference names and lengths.

    @param filenameOrSamfile: Either a C{str} SAM/BAM file name or an
        instance of C{pysam.AlignmentFile}.
    @param indent: An C{int} number of spaces to indent each line.
    @return: A C{str} describing known reference names and their lengths.
    """
    indent = ' ' * indent

    def _references(sam):
        result = []
        for i in range(sam.nreferences):
            result.append('%s%s (length %d)' % (
                indent, sam.get_reference_name(i), sam.lengths[i]))
        return '\n'.join(result)

    if isinstance(filenameOrSamfile, six.string_types):
        with samfile(filenameOrSamfile) as sam:
            return _references(sam)
    else:
        return _references(sam)


class SAMFilter(object):
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
    def __init__(self, filename, filterRead=None, referenceIds=None,
                 storeQueryIds=False, dropUnmapped=False,
                 dropSecondary=False, dropSupplementary=False,
                 dropDuplicates=False, keepQCFailures=False, minScore=None,
                 maxScore=None, scoreTag='AS'):
        self.filename = filename
        self.filterRead = filterRead
        self.referenceIds = referenceIds
        self.dropUnmapped = dropUnmapped
        self.dropSecondary = dropSecondary
        self.dropSupplementary = dropSupplementary
        self.dropDuplicates = dropDuplicates
        self.keepQCFailures = keepQCFailures
        self.storeQueryIds = storeQueryIds
        self.minScore = minScore
        self.maxScore = maxScore
        self.scoreTag = scoreTag

    @staticmethod
    def addFilteringOptions(parser, samfileIsPositionalArg=False):
        """
        Add options to an argument parser for filtering SAM/BAM.

        @param samfileIsPositionalArg: If C{True} the SAM/BAM file must
            be given as the final argument on the command line (without
            being preceded by --samfile).
        @param parser: An C{argparse.ArgumentParser} instance.
        """
        parser.add_argument(
            '%ssamfile' % ('' if samfileIsPositionalArg else '--'),
            required=True,
            help='The SAM/BAM file to filter.')

        parser.add_argument(
            '--referenceId', metavar='ID', nargs='+', action='append',
            help=('A reference sequence id whose alignments should be kept '
                  '(alignments against other references will be dropped). '
                  'If omitted, alignments against all references will be '
                  'kept. May be repeated.'))

        parser.add_argument(
            '--dropUnmapped', default=False, action='store_true',
            help='If given, unmapped matches will not be output.')

        parser.add_argument(
            '--dropSecondary', default=False, action='store_true',
            help='If given, secondary matches will not be output.')

        parser.add_argument(
            '--dropSupplementary', default=False, action='store_true',
            help='If given, supplementary matches will not be output.')

        parser.add_argument(
            '--dropDuplicates', default=False, action='store_true',
            help=('If given, matches flagged as optical or PCR duplicates '
                  'will not be output.'))

        parser.add_argument(
            '--keepQCFailures', default=False, action='store_true',
            help=('If given, reads that are considered quality control '
                  'failures will be included in the output.'))

        parser.add_argument(
            '--minScore', type=float,
            help=('If given, alignments with --scoreTag (default AS) values '
                  'less than this value will not be output. If given, '
                  'alignments that do not have a score will not be output.'))

        parser.add_argument(
            '--maxScore', type=float,
            help=('If given, alignments with --scoreTag (default AS) values '
                  'greater than this value will not be output. If given, '
                  'alignments that do not have a score will not be output.'))

        parser.add_argument(
            '--scoreTag', default='AS',
            help=('The alignment tag to extract for --minScore and --maxScore '
                  'comparisons.'))

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
        referenceIds = (set(chain.from_iterable(args.referenceId))
                        if args.referenceId else None)

        return cls(
            args.samfile,
            filterRead=filterRead,
            referenceIds=referenceIds,
            storeQueryIds=storeQueryIds,
            dropUnmapped=args.dropUnmapped,
            dropSecondary=args.dropSecondary,
            dropSupplementary=args.dropSupplementary,
            dropDuplicates=args.dropDuplicates,
            keepQCFailures=args.keepQCFailures,
            minScore=args.minScore,
            maxScore=args.maxScore)

    def alignments(self):
        """
        Get alignments from the SAM/BAM file, subject to filtering.
        """
        referenceIds = self.referenceIds
        dropUnmapped = self.dropUnmapped
        dropSecondary = self.dropSecondary
        dropSupplementary = self.dropSupplementary
        dropDuplicates = self.dropDuplicates
        keepQCFailures = self.keepQCFailures
        storeQueryIds = self.storeQueryIds
        filterRead = self.filterRead
        minScore = self.minScore
        maxScore = self.maxScore
        scoreTag = self.scoreTag

        if storeQueryIds:
            self.queryIds = queryIds = set()

        count = 0
        with samfile(self.filename) as samAlignment:
            for count, alignment in enumerate(samAlignment.fetch(), start=1):
                if storeQueryIds:
                    queryIds.add(alignment.query_name)

                if minScore is not None or maxScore is not None:
                    try:
                        score = alignment.get_tag(scoreTag)
                    except KeyError:
                        continue
                    else:
                        if ((minScore is not None and score < minScore) or
                                (maxScore is not None and score > maxScore)):
                            continue

                if ((filterRead is None or
                     filterRead(Read(alignment.query_name,
                                     alignment.query_sequence,
                                     alignment.qual))) and
                    not (
                        (referenceIds and
                         alignment.reference_name not in referenceIds) or
                        (alignment.is_unmapped and dropUnmapped) or
                        (alignment.is_secondary and dropSecondary) or
                        (alignment.is_supplementary and dropSupplementary) or
                        (alignment.is_duplicate and dropDuplicates) or
                        (alignment.is_qcfail and not keepQCFailures))):
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
                            'Reference %r is not present in the SAM/BAM file.'
                            % referenceId)
                    else:
                        result[referenceId] = sam.lengths[tid]
            else:
                result = dict(zip(sam.references, sam.lengths))

        return result


class PaddedSAM(object):
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
                'Your %d SAM/BAM file reference sequence '
                'lengths (%s) are not all identical.' % (
                    len(referenceLengths),
                    ', '.join(
                        '%s=%d' % (id_, referenceLengths[id_])
                        for id_ in sorted(referenceLengths))))

        # Get the length of any of the sequences (they are all identical).
        self.referenceLength = referenceLengths.popitem()[1]
        self.samFilter = samFilter
        # self.referenceInsertions will be keyed by query id (the query
        # that would cause a reference insertion). The values will be lists
        # of 2-tuples, with each 2-tuple containing an offset into the
        # reference sequence and the C{str} of nucleotides that would be
        # inserted starting at that offset.
        self.referenceInsertions = defaultdict(list)

    def queries(self, rcSuffix='', rcNeeded=False, padChar='-',
                queryInsertionChar='N', allowDuplicateIds=False,
                addAlignment=False):
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
        lastQuery = None

        for lineNumber, alignment in enumerate(
                self.samFilter.alignments(), start=1):

            query = alignment.query_sequence

            # Secondary (and presumably supplementary) alignments may have
            # a '*' (None in pysam) SEQ field, indicating that the previous
            # sequence should be used. This is best practice according to
            # section 2.5.2 of https://samtools.github.io/hts-specs/SAMv1.pdf
            if query is None:
                if alignment.is_secondary or alignment.is_supplementary:
                    if lastQuery is None:
                        raise InvalidSAM(
                            'Query line %d has an empty SEQ field, but no '
                            'previous alignment is present.' % lineNumber)
                    else:
                        query = lastQuery
                else:
                    raise InvalidSAM(
                        'Query line %d has an empty SEQ field, but the '
                        'alignment is not marked as secondary or '
                        'supplementary.' % lineNumber)

            # Remember the last query here (before we potentially modify it
            # due to it being reverse complimented for the alignment).
            lastQuery = query

            if alignment.is_reverse:
                if rcNeeded:
                    query = DNARead('id', query).reverseComplement().sequence
                if rcSuffix:
                    alignment.query_name += rcSuffix

            # Adjust the query id if it's a duplicate and we're not allowing
            # duplicates.
            if allowDuplicateIds:
                queryId = alignment.query_name
            else:
                count = idCount[alignment.query_name]
                idCount[alignment.query_name] += 1
                queryId = alignment.query_name + (
                    '' if count == 0 else '/%d' % count)

            referenceStart = alignment.reference_start
            atStart = True
            queryIndex = 0
            referenceIndex = referenceStart
            alignedSequence = ''

            for operation, length in alignment.cigartuples:

                # The operations are tested in the order they appear in
                # https://samtools.github.io/hts-specs/SAMv1.pdf It would be
                # more efficient to test them in order of frequency of
                # occurrence.
                if operation in MATCH_OPERATIONS:
                    atStart = False
                    alignedSequence += query[queryIndex:queryIndex + length]
                elif operation == CINS:
                    # Insertion to the reference. This consumes query bases but
                    # we don't output them because the reference cannot be
                    # changed.  I.e., these bases in the query would need to be
                    # inserted into the reference.  Remove these bases from the
                    # query but record what would have been inserted into the
                    # reference.
                    atStart = False
                    self.referenceInsertions[queryId].append(
                        (referenceIndex,
                         query[queryIndex:queryIndex + length]))
                elif operation == CDEL:
                    # Delete from the reference. Some bases from the reference
                    # would need to be deleted to continue the match. So we put
                    # an insertion into the query to compensate.
                    atStart = False
                    alignedSequence += queryInsertionChar * length
                elif operation == CREF_SKIP:
                    # Skipped reference. Opens a gap in the query. For
                    # mRNA-to-genome alignment, an N operation represents an
                    # intron.  For other types of alignments, the
                    # interpretation of N is not defined. So this is unlikely
                    # to occur.
                    atStart = False
                    alignedSequence += queryInsertionChar * length
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
                            alignedSequence += query[queryIndex + unwantedLeft:
                                                     queryIndex + length]
                            referenceStart = 0
                        else:
                            referenceStart -= length
                            alignedSequence += query[
                                queryIndex:queryIndex + length]
                    else:
                        unwantedRight = (
                            (referenceStart + len(alignedSequence) + length) -
                            referenceLength)

                        if unwantedRight > 0:
                            # The query protrudes right. Copy its left part.
                            alignedSequence += query[
                                queryIndex:queryIndex + length - unwantedRight]
                        else:
                            alignedSequence += query[
                                queryIndex:queryIndex + length]
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
                    raise ValueError('Unknown CIGAR operation:', operation)

                if operation in _CONSUMES_QUERY:
                    queryIndex += length

                if operation in _CONSUMES_REFERENCE:
                    referenceIndex += length

            # Sanity check that we consumed the entire query.
            assert queryIndex == len(query)

            # We cannot test we consumed the entire reference.  The CIGAR
            # string applies to (and exhausts) the query but is silent
            # about the part of the reference that lies to the right of the
            # aligned query.

            # Put gap characters before and after the aligned sequence so that
            # it is offset properly and matches the length of the reference.
            paddedSequence = (
                (padChar * referenceStart) +
                alignedSequence +
                padChar * (referenceLength -
                           (referenceStart + len(alignedSequence))))

            read = Read(queryId, paddedSequence)
            if addAlignment:
                read.alignment = alignment

            yield read
