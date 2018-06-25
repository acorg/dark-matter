from contextlib import contextmanager
from collections import Counter, defaultdict

from pysam import (
    AlignmentFile, CMATCH, CINS, CDEL, CREF_SKIP, CSOFT_CLIP, CHARD_CLIP, CPAD,
    CEQUAL, CDIFF)

from dark.reads import Read, DNARead


class UnequalReferenceLengthError(Exception):
    "The reference sequences in a SAM/BAM file are not all of the same length."


class UnknownReference(Exception):
    "Reference sequence not found in SAM/BAM file."


# From https://samtools.github.io/hts-specs/SAMv1.pdf
_CONSUMES_QUERY = {CMATCH, CINS, CSOFT_CLIP, CEQUAL, CDIFF}
_CONSUMES_REFERENCE = {CMATCH, CDEL, CREF_SKIP, CEQUAL, CDIFF}


class PaddedSAM(object):
    """
    Obtain aligned (padded) queries from a SAM/BAM file.

    @param filename: The C{str} name of the SAM/BAM file.
    """
    def __init__(self, filename):
        self.samfile = AlignmentFile(filename)
        # self.referenceInsertions will be keyed by query id (the query
        # that would cause a reference insertion). The values will be lists
        # of 2-tuples, with each 2-tuple containing an offset into the
        # reference sequence and the C{str} of nucleotide that would be
        # inserted starting at that offset.
        self.referenceInsertions = defaultdict(list)

    def close(self):
        """
        Close the opened SAM/BAM file.
        """
        self.samfile.close()

    def referencesToStr(self, indent=0):
        """
        List the reference names and their lengths.

        @param indent: An C{int} number of spaces to indent each line.
        @return: A C{str} describing known reference names and their lengths.
        """
        samfile = self.samfile
        result = []
        indent = ' ' * indent
        for i in range(samfile.nreferences):
            result.append('%s%s (length %d)' % (
                indent, samfile.get_reference_name(i), samfile.lengths[i]))
        return '\n'.join(result)

    def queries(self, referenceName=None, minLength=0, rcSuffix='',
                dropSecondary=False, dropSupplementary=False,
                dropDuplicates=False, allowDuplicateIds=False,
                keepQCFailures=False, rcNeeded=False, padChar='-',
                queryInsertionChar='N'):
        """
        Produce padded (with gaps) queries according to the CIGAR string and
        reference sequence length for each matching query sequence.

        @param referenceName: The C{str} name of the reference sequence to
            print alignments for. This is only needed if the SAM/BAM alignment
            was made against multiple references *and* they are not all of the
            same length. If there is only one reference sequence or if all
            reference sequences are of the same length, there is no need to
            provide a reference name (i.e., pass C{None}).
        @param minLength: Ignore queries shorter than this C{int} value. Note
            that this refers to the length of the query sequence once it has
            been aligned to the reference. The alignment may introduce
            C{queryInsertionChar} characters into the read, and these are
            counted towards its length because the alignment is assuming the
            query is missing a base at those locations.
        @param rcSuffix: A C{str} to add to the end of query names that are
            reverse complemented. This is added before the /1, /2, etc., that
            are added for duplicated ids (if there are duplicates and
            C{allowDuplicateIds} is C{False}.
        @param dropSecondary: If C{True}, secondary matches will not be
            yielded.
        @param dropSupplementary: If C{True}, supplementary matches will not be
            yielded.
        @param dropDuplicates: If C{True}, matches flagged as optical or PCR
            duplicates will not be yielded.
        @param allowDuplicateIds: If C{True}, repeated query ids (due to
            secondary or supplemental matches) will not have /1, /2, etc.
            appended to their ids. So repeated ids may appear in the yielded
            FASTA.
        @param keepQCFailures: If C{True}, reads that are marked as quality
            control failures will be included in the output.
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
        @raises UnequalReferenceLengthError: If C{referenceName} is C{None}
            and the reference sequence lengths in the SAM/BAM file are not all
            identical.
        @raises UnknownReference: If C{referenceName} does not exist.
        @return: A generator that yields C{Read} instances that are padded
            with gap characters to align them to the length of the reference
            sequence.
        """
        samfile = self.samfile

        if referenceName:
            referenceId = samfile.get_tid(referenceName)
            if referenceId == -1:
                raise UnknownReference(
                    'Reference %r is not present in the SAM/BAM file.'
                    % referenceName)
            referenceLength = samfile.lengths[referenceId]
        else:
            # No reference given. All references must have the same length.
            if len(set(samfile.lengths)) != 1:
                raise UnequalReferenceLengthError(
                    'Your SAM/BAM file has %d reference sequences, and their '
                    'lengths (%s) are not all identical.' % (
                        samfile.nreferences,
                        ', '.join(map(str, sorted(samfile.lengths)))))
            referenceId = None
            referenceLength = samfile.lengths[0]

        # Hold the count for each id so we can add /1, /2 etc to duplicate
        # ids (unless --allowDuplicateIds was given).
        idCount = Counter()

        MATCH_OPERATIONS = {CMATCH, CEQUAL, CDIFF}
        lastQuery = None

        for lineNumber, read in enumerate(samfile.fetch(), start=1):
            if (read.is_unmapped or
                    (read.is_secondary and dropSecondary) or
                    (read.is_supplementary and dropSupplementary) or
                    (read.is_duplicate and dropDuplicates) or
                    (read.is_qcfail and not keepQCFailures) or
                    (referenceId is not None and
                     read.reference_id != referenceId)):
                continue

            query = read.query_sequence

            # Secondary (and presumably supplementary) alignments may have
            # a '*' (None in pysam) SEQ field, indicating that the previous
            # sequence should be used. This is best practice according to
            # section 2.5.2 of https://samtools.github.io/hts-specs/SAMv1.pdf
            if query is None:
                if read.is_secondary or read.is_supplementary:
                    if lastQuery is None:
                        raise ValueError(
                            'Query line %d has an empty SEQ field, but no '
                            'previous alignment is present.' % lineNumber)
                    else:
                        query = lastQuery
                else:
                    raise ValueError(
                        'Query line %d has an empty SEQ field, but the '
                        'alignment is not marked as secondary or '
                        'supplementary.' % lineNumber)

            # Remember the last query here (before we potentially modify it
            # due to it being reverse complimented for the alignment).
            lastQuery = query

            if read.is_reverse:
                if rcNeeded:
                    query = DNARead('id', query).reverseComplement().sequence
                if rcSuffix:
                    read.query_name += rcSuffix

            # Adjust the query id if it's a duplicate and we're not
            # allowing duplicates.
            if allowDuplicateIds:
                queryId = read.query_name
            else:
                count = idCount[read.query_name]
                idCount[read.query_name] += 1
                queryId = read.query_name + (
                    '' if count == 0 else '/%d' % count)

            referenceStart = read.reference_start
            atStart = True
            queryIndex = 0
            referenceIndex = referenceStart
            alignedSequence = ''

            for operation, length in read.cigartuples:

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
            # string applies to (i.e., exhausts) the query and is silent about
            # the part of the reference that to the right of the aligned query.

            # Check the length restriction now that we have (possibly) added
            # queryInsertionChar characters to pad the query out to the length
            # it requires to match the reference.
            if len(alignedSequence) < minLength:
                continue

            # Put gap characters before and after the aligned sequence so that
            # it is offset properly and matches the length of the reference.
            paddedSequence = (
                (padChar * referenceStart) +
                alignedSequence +
                padChar * (referenceLength -
                           (referenceStart + len(alignedSequence))))

            yield Read(queryId, paddedSequence)


@contextmanager
def samfile(filename):
    """
    A context manager to open and close a SAM/BAM file.

    @param filename: A C{str} file name to open.
    """
    f = AlignmentFile(filename)
    yield f
    f.close()
