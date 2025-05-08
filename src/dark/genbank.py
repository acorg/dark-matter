from warnings import warn


class GenomeRanges:
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
        if rangeStr.startswith("join{") and rangeStr[-1] == "}":
            join = True
            inner = rangeStr[5:-1]
        else:
            join = False
            inner = rangeStr

        ranges = []
        subRanges = inner.split(", ")
        nRanges = len(subRanges)

        if nRanges == 1 and join:
            raise ValueError(
                'Could not parse GenBank range string "%s". '
                "join{} can only be used with multiple ranges." % rangeStr
            )
        elif nRanges > 1 and not join:
            raise ValueError(
                'Could not parse GenBank range string "%s". '
                "Multiple ranges must be wrapped in join{}." % rangeStr
            )

        for subRange in subRanges:
            if subRange.endswith("](+)"):
                forward = True
            elif subRange.endswith("](-)"):
                forward = False
            else:
                raise ValueError(
                    'Could not parse GenBank range string "%s". '
                    'Range "%s" does not end with ](+) or ](-).' % (rangeStr, subRange)
                )

            if not subRange.startswith("["):
                raise ValueError(
                    'Could not parse GenBank range string "%s". '
                    'Range "%s" does not start with "[".' % (rangeStr, subRange)
                )
            try:
                start, stop = subRange[1:-4].split(":")
            except ValueError as e:
                raise ValueError(
                    'Could not parse GenBank range string "%s". '
                    'Original parsing ValueError was "%s".' % (rangeStr, e)
                )
            else:
                if start.startswith("<"):
                    start = start[1:]
                if stop.startswith(">"):
                    stop = stop[1:]
                start, stop = map(int, (start, stop))
                if start > stop:
                    raise ValueError(
                        'Could not parse GenBank range string "%s". '
                        "Offset values (%d, %d) cannot decrease."
                        % (rangeStr, start, stop)
                    )

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
                    warn(
                        "Contiguous GenBank ranges detected: [%d:%d] "
                        "followed by [%d:%d]." % (lastStart, lastStop, start, stop)
                    )
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
        return "<%s: %s>" % (self.__class__.__name__, ", ".join(map(str, self.ranges)))

    def circular(self, genomeLength):
        """
        Determine whether the offset ranges of a protein in a genome span the
        end of the genome (indicating that the genome may be circular).

        @param genomeLength: The C{int} length of the genome.
        @return: A C{bool}, C{True} if the ranges overlap the genome end.
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
            # ((0, 100, False),) to indicate circularity is a bit arbitrary.
            # Here I've decided that it should not because I intend to use
            # this function to decide whether to output multiple SAM lines
            # for the match of a read against a reference. In this
            # degenerate case just a single SAM line suffices, whereas in a
            # non-degenerate case (such as ((0, 75, False), (75, 100, False))
            # on a genome of length 100), a chimeric SAM alignment must be
            # used, which requires quite different processing than a normal
            # linear match. YMMV, and if you don't like it you can easily
            # add an argument to this function that changes this behaviour
            # (leaving the default to do what the code currently does).
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
        offsetInGenome = remaining = (match["sstart"] - 1) * 3

        for start, stop, _ in self.ranges:
            rangeWidth = stop - start
            if remaining < rangeWidth:
                # The match starts in this range.
                return start + remaining
            else:
                remaining -= rangeWidth
        else:
            raise ValueError(
                "Starting nucleotide offset %d not found in protein "
                "nucleotide ranges %s."
                % (
                    offsetInGenome,
                    ", ".join(("(%d, %d)" % (i, j)) for i, j, _ in self.ranges),
                )
            )

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


def getSourceInfo(genome, keys=("host", "note", "organism", "mol_type"), logfp=None):
    """
    Extract summary information from a genome source feature.

    @param genome: A GenBank genome record, as parsed by SeqIO.parse
    @param keys: An iterable of C{str} keys to extract from the Genbank
        source features.
    @param logfp: If not C{None}, a file pointer to write verbose
        progress output to.
    @return: A C{dict} with keys for the various pieces of information
        (if any) found in the source feature (see the return value below
        for detail). Or C{None} if no source feature is found or a source
        feature does not have length 1.
    """
    result = {}

    for feature in genome.features:
        if feature.type == "source":
            for key in keys:
                try:
                    values = feature.qualifiers[key]
                except KeyError:
                    value = None
                    if key != "note" and logfp:
                        print(
                            "Genome %r (accession %s) source info has "
                            "no %r feature." % (genome.description, genome.id, key),
                            file=logfp,
                        )
                else:
                    if len(values) == 1:
                        value = values[0]

                        if key == "mol_type":
                            assert value[-3:] in ("DNA", "RNA")

                    elif len(values) > 1 and key == "host":
                        value = ", ".join(values)
                    else:
                        if logfp:
                            print(
                                "Genome %r (accession %s) has source "
                                "feature %r with length != 1: %r"
                                % (genome.description, genome.id, key, values),
                                file=logfp,
                            )
                        return

                result[key] = value

            # Break so no other features are examined; we only want source.
            break
    else:
        if logfp:
            print(
                "Genome %r (accession %s) had no source feature! "
                "Skipping." % (genome.description, genome.id),
                file=logfp,
            )
        return

    return result


def getCDSInfo(
    genome, feature, keys=("gene", "note", "product", "protein_id", "translation")
):
    """
    Extract summary information from a genome CDS feature.

    @param genome: A GenBank genome record, as parsed by SeqIO.parse
    @param feature: A feature from a genome, as produced by BioPython's
        GenBank parser.
    @param keys: An iterable of C{str} keys to extract from the Genbank
        feature.
    @return: A C{dict} with keys for the various pieces of information
        found in the feature (see the return value below for detail).
        Or C{None} if the feature is not of interest or otherwise invalid.
    """
    qualifiers = feature.qualifiers

    # Check in advance that all feature qualifiers we're interested in
    # have the right lengths, if they're present.
    for key in "gene", "note", "product", "protein_id", "translation":
        if key in qualifiers:
            assert (
                len(qualifiers[key]) == 1
            ), "GenBank qualifier key %s is not length one %r" % (key, qualifiers[key])

    # A protein id is mandatory.
    if "protein_id" in qualifiers:
        proteinId = qualifiers["protein_id"][0]
    else:
        if "translation" in qualifiers:
            warn(
                "Genome %r (accession %s) has CDS feature with no "
                "protein_id feature but has a translation! "
                "Skipping.\nFeature: %s" % (genome.description, genome.id, feature)
            )
        return

    # A translated (i.e., amino acid) sequence is mandatory.
    if "translation" in qualifiers:
        translation = qualifiers["translation"][0]
    else:
        warn(
            "Genome %r (accession %s) has CDS feature with protein "
            "%r with no translated sequence. Skipping."
            % (genome.description, genome.id, proteinId)
        )
        return

    featureLocation = str(feature.location)

    # Make sure the feature's location string can be parsed.
    try:
        ranges = GenomeRanges(featureLocation)
    except ValueError as e:
        warn(
            "Genome %r  (accession %s) contains unparseable CDS "
            "location for protein %r. Skipping. Error: %s"
            % (genome.description, genome.id, proteinId, e)
        )
        return
    else:
        # Does the protein span the end of the genome? This indicates a
        # circular genome.
        circular = int(ranges.circular(len(genome.seq)))

    if feature.location.start >= feature.location.end:
        warn(
            "Genome %r (accession %s) contains feature with start "
            "(%d) >= stop (%d). Skipping.\nFeature: %s"
            % (
                genome.description,
                genome.id,
                feature.location.start,
                feature.location.end,
                feature,
            )
        )
        return

    strand = feature.strand
    if strand is None:
        # The strands of the protein in the genome are not all the same
        # (see Bio.SeqFeature.CompoundLocation._get_strand).  The
        # protein is formed by the combination of reading one strand in
        # one direction and the other in the other direction.
        #
        # This occurs just once in all 1.17M proteins found in all 700K
        # RVDB (C-RVDBv15.1) genomes, for protein YP_656697.1 on the
        # Ranid herpesvirus 1 strain McKinnell genome (NC_008211.1).
        #
        # This situation makes turning DIAMOND protein output into
        # SAM very complicated because a match on such a protein
        # cannot be stored as a SAM linear alignment. It instead
        # requires a multi-line 'supplementary' alignment. The code
        # and tests for that are more complex than I want to deal
        # with at the moment, just for the sake of one protein in a
        # frog herpesvirus.
        warn(
            "Genome %s (accession %s) has protein %r with mixed "
            "orientation!" % (genome.description, genome.id, proteinId)
        )
        return
    elif strand == 0:
        # This never occurs for proteins corresponding to genomes in
        # the RVDB database C-RVDBv15.1.
        warn(
            "Genome %r (accession %s) has protein %r with feature "
            "with strand of zero!" % (genome.description, genome.id, proteinId)
        )
        return
    else:
        assert strand in (1, -1)
        forward = strand == 1
        # Make sure the strand agrees with the orientations in the
        # string BioPython makes out of the locations.
        assert ranges.orientations() == {forward}

    return {
        "circular": circular,
        "featureLocation": featureLocation,
        "forward": forward,
        "gene": qualifiers.get("gene", [""])[0],
        "note": qualifiers.get("note", [""])[0],
        "product": qualifiers.get("product", ["UNKNOWN"])[0],
        "proteinId": proteinId,
        "ranges": ranges,
        "strand": strand,
        "translation": translation,
    }
