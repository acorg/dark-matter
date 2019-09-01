import warnings


class GenomeRanges(object):
    """
    Split and manage a GenBank range string as converted by BioPython.

    @param rangeStr: A C{str} indicating the (0-based) genome nucleotide
        offsets covered by a protein. This can have the following example
        forms:

        [9462:10137](+)
        [11969:12575](-)
        join{[126386:126881](-), [125941:126232](+)}
        join{[153005:153821](+), [82030:82618](-), [80480:81305](-)}

        Note that this is the string as returned from BioPython SeqIO.parse
        when it parses a GenBank flat file, not the string that is in such a
        file. For that use the C{splitRange} function you can find in an
        earlier version of this file (i.e., look in its git history).

    @raise ValueError: If C{rangeStr} cannot be correctly parsed.
    @return: A C{list} whose elements are 3-C{tuple}s holding the C{int} start
        and stop offsets and a C{bool} that indicates whether the offsets
        corresponds to the forward strand (i.e., will be C{True} when
        there is a '(+)' in a range and C{False} when there is a '(-)' which
        indicates the reverse strand).

        Example arguments and their return values:

            '[9462:10137](+)' => ((9462, 10137, True),)

            '[11969:12575](-)' => ((11969, 12575, False),)

            'join{[126386:126881](-), [125941:126232](+)}' =>
                ((126386, 126881, False),
                 (125941, 126232, True))
    """
    def __init__(self, rangeStr):
        if rangeStr.startswith('join{') and rangeStr[-1] == '}':
            join = True
            inner = rangeStr[5:-1]
        else:
            join = False
            inner = rangeStr

        ranges = []
        subRanges = inner.split(', ')
        nRanges = len(subRanges)

        if nRanges == 1 and join:
            raise ValueError('Could not parse GenBank range string "%s". '
                             'join{} can only be used with multiple ranges.' %
                             rangeStr)
        elif nRanges > 1 and not join:
            raise ValueError('Could not parse GenBank range string "%s". '
                             'Multiple ranges must be wrapped in join{}.' %
                             rangeStr)

        for subRange in subRanges:
            if subRange.endswith('](+)'):
                forward = True
            elif subRange.endswith('](-)'):
                forward = False
            else:
                raise ValueError('Could not parse GenBank range string "%s". '
                                 'Range "%s" does not end with ](+) or ](-).' %
                                 (rangeStr, subRange))

            if not subRange.startswith('['):
                raise ValueError('Could not parse GenBank range string "%s". '
                                 'Range "%s" does not start with "[".' %
                                 (rangeStr, subRange))
            try:
                start, stop = subRange[1:-4].split(':')
            except ValueError as e:
                raise ValueError('Could not parse GenBank range string "%s". '
                                 'Original parsing ValueError was "%s".' %
                                 (rangeStr, e))
            else:
                if start.startswith('<'):
                    start = start[1:]
                if stop.startswith('>'):
                    stop = stop[1:]
                start, stop = map(int, (start, stop))
                if start > stop:
                    raise ValueError(
                        'Could not parse GenBank range string "%s". '
                        'Offset values (%d, %d) cannot decrease.' %
                        (rangeStr, start, stop))

                ranges.append((start, stop, forward))

        self.ranges = self._mergeContiguousRanges(ranges)
        self._nRanges = len(self.ranges)

    def _mergeContiguousRanges(self, ranges):
        """
        Merge ranges that are contiguous (follow each other immediately on the
        genome).

        @param ranges: An iterable of (start, stop, forward) tuples.
        @return: A C{tuple} of (start, stop, forward) tuples with contiguous
            ranges found in C{ranges} merged.
        """
        result = []
        lastStart = lastStop = lastForward = None

        for index, (start, stop, forward) in enumerate(ranges):
            if lastStart is None:
                lastStart, lastStop, lastForward = start, stop, forward
            else:
                if start == lastStop and forward == lastForward:
                    # This range continues the previous one.
                    warnings.warn(
                        'Contiguous GenBank ranges detected: [%d:%d] '
                        'followed by [%d:%d].' %
                        (lastStart, lastStop, start, stop))
                    lastStop = stop
                else:
                    # Emit the range that just got terminated.
                    result.append((lastStart, lastStop, lastForward))
                    # And remember this one.
                    lastStart, lastStop, lastForward = start, stop, forward

        # Emit the final range.
        result.append((lastStart, lastStop, lastForward))

        return tuple(result)

    def __str__(self):
        return '<%s: %s>' % (
            self.__class__.__name__, ', '.join(map(str, self.ranges)))

    def circular(self, genomeLength):
        """
        Determine whether the offset ranges of a protein in a genome span the
        end of the genome (indicating that the genome may be circular).

        @param genomeLength: The C{int} length of the genome.
        @return: A C{bool}, C{True} if the ranges overlaps the genome end.
        """
        if self._nRanges == 1:
            # If there is only one range, we simply return False even though it
            # could be that the single range spans the whole genome. That is
            # not the same as having two ranges, one of which ends at the end
            # of the genome and the next starting at the start of the genome.
            #
            # I.e., a range tuples like this ((0, 100, False),) on a genome of
            # length 100 will return False, whereas ((0, 75, False), (75, 100,
            # False)) on a genome of length 100 will return True.
            #
            # The decision about whether you consider the degenerate case of
            # ((0, 100, False),) to indicate circularity is a bit
            # arbitrary. Here I've decided that it should not because I intend
            # to use this function to decide whether to output multiple SAM
            # lines for the match of a read against a reference. In this
            # degenerate case just a single SAM line suffices, whereas in a
            # non-degenerate case (such as ((0, 75, False), (75, 100, False))
            # on a genome of length 100), a chimeric SAM alignment must be
            # used, which requires quite different processing than a normal
            # linear match. YMMV, and if you don't like it you can easily add
            # an argument to this function that changes this behaviour (leaving
            # the default to do what the code currently does).
            return False

        for index, (start, stop, _) in enumerate(self.ranges):
            if stop == genomeLength:
                return self.ranges[(index + 1) % self._nRanges][0] == 0

        return False

    def startInGenome(self, match):
        """
        Calculate the position in the nucleotide genome where a DIAMOND match
        begins.

        @param match: A C{dict} with information about the DIAMOND match, as
            returned by C{self._preprocessMatch}. Must contain a 'sstart' key
            giving a 1-based C{int} offset of the start of the match in the
            protein.
        @return: The C{int} offset of the start of the DIAMOND match in the
            genome.
        """

        # Calculate the start of the match in the genome, given its start in
        # the protein.
        offsetInGenome = remaining = (match['sstart'] - 1) * 3

        for start, stop, _ in self.ranges:
            rangeWidth = stop - start
            if remaining < rangeWidth:
                # The match starts in this range.
                return start + remaining
            else:
                remaining -= rangeWidth
        else:
            raise ValueError(
                'Starting nucleotide offset %d not found in protein '
                'nucleotide ranges %s.' %
                (offsetInGenome,
                 ', '.join(('(%d, %d)' % (i, j))
                           for i, j, _ in self.ranges)))

    def orientations(self):
        """
        Produce the set of all orientations for our ranges.

        @return: A C{set} of C{True} and C{False} range orientations.
        """
        return set(r[2] for r in self.ranges)

    def distinctRangeCount(self, genomeLength):
        """
        Determine the number of distinct ranges, given the genome length.

        @param genomeLength: The C{int} length of the genome.
        @return: An C{int} which is the number of ranges that are not
            contiguous with one another.
        """
        return len(self.ranges) - int(self.circular(genomeLength))
