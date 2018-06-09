from collections import Counter

from pysam import (
    AlignmentFile, CMATCH, CINS, CDEL, CREF_SKIP, CSOFT_CLIP, CHARD_CLIP, CPAD,
    CEQUAL, CDIFF)

from dark.reads import Read, DNARead


class UnequalReferenceLengthError(Exception):
    "The reference sequences in a SAM/BAM file are not all of the same length."


class UnknownReference(Exception):
    "Reference sequence not found in SAM/BAM file."


class PaddedSAM(object):
    """
    Obtain aligned (padded) queries from a SAM/BAM file.

    @param filename: The C{str} name of the SAM/BAM file.
    """
    def __init__(self, filename):
        self.samfile = AlignmentFile(filename)

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
                allowDuplicateIds=False, keepQCFailures=False, rcNeeded=False):
        """
        Produce padded (with '-') queries according to the CIGAR string and
        reference sequence length for each matching query sequence.

        @param referenceName: The C{str} name of the reference sequence to
            print alignments for. This is only needed if the SAM/BAM alignment
            was made against multiple references *and* they are not all of the
            same length. If there is only one reference sequence or if all
            reference sequences are of the same length, there is no need to
            provide a reference name (i.e., pass C{None}).
        @param minLength: Ignore queries shorter than this C{int} value.
        @param rcSuffix: A C{str} to add to the end of query names that are
            reverse complemented. This is added before the /1, /2, etc., that
            are added for duplicated ids (if there are duplicates and
            C{allowDuplicateIds} is C{False}.
        @param dropSecondary: If C{True}, secondary matches will not be
            yielded.
        @param dropSupplementary: If C{True}, supplementary matches will not be
            yielded.
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
        @raises UnequalReferenceLengthError: If C{referenceName} is C{None}
            and the reference sequence lengths in the SAM/BAM file are not all
            identical.
        @raises UnknownReference: If C{referenceName} does not exist.
        @return: A generator that yields C{Read} instances that are padded
            with '-' characters to align them to the length of the reference
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

        for read in samfile.fetch():
            query = read.query_sequence
            queryLength = len(query)
            if (read.is_unmapped or
                    queryLength < minLength or
                    (read.is_secondary and dropSecondary) or
                    (read.is_supplementary and dropSupplementary) or
                    (read.is_qcfail and not keepQCFailures) or
                    (referenceId is not None and
                     read.reference_id != referenceId)):
                continue

            if rcNeeded and read.is_reverse:
                query = DNARead('id', query).reverseComplement().sequence
                if rcSuffix:
                    read.query_name += rcSuffix

            index = 0
            alignedSequence = []
            for operation, length in read.cigartuples:

                if operation in MATCH_OPERATIONS:
                    alignedSequence.append(query[index:index + length])
                    index += length
                elif operation == CINS:
                    # Insertion to the reference. This consumes query bases
                    # but we don't output them because the reference cannot
                    # be changed.  I.e., these bases in the query would
                    # need to be inserted into the reference.
                    #
                    # But what do we do? Removing query bases from the
                    # middle of the query seems wrong.
                    index += length
                elif operation == CDEL:
                    # Delete from the reference. Some bases from the
                    # reference would need to be deleted to continue the
                    # match. So we put gaps into the query.
                    alignedSequence.append('-' * length)
                elif operation == CREF_SKIP:
                    # Is this the same as CDEL (in terms of what we do)?
                    pass
                elif operation == CSOFT_CLIP:
                    # Bases in the query that are not part of the match. We
                    # remove these from the query. I think this means these
                    # bases protrude, either before the start or after the
                    # end of the reference.
                    #
                    # S may only have H operations between them and the
                    # ends of the CIGAR string.
                    index += length
                elif operation == CHARD_CLIP:
                    # Some bases have been completely removed from the
                    # query. Can only be present as the first and/or last
                    # operation.
                    pass
                elif operation == CPAD:
                    pass
                else:
                    raise ValueError('Unknown CIGAR operation:', operation)

            # Sanity check that we consumed the entire query.
            assert index == queryLength

            alignedSequence = ''.join(alignedSequence)

            # Put '-' gap characters before and after the aligned sequence
            # so that it is offset properly and matches the length of the
            # reference.
            paddedSequence = (
                ('-' * read.reference_start) +
                alignedSequence +
                '-' * (referenceLength -
                       (read.reference_start + len(alignedSequence))))

            if allowDuplicateIds:
                suffix = ''
            else:
                count = idCount[read.query_name]
                idCount[read.query_name] += 1
                suffix = '' if count == 0 else '/%d' % count

            yield Read('%s%s' % (read.query_name, suffix), paddedSequence)
