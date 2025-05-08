from Bio.Seq import reverse_complement


def findPrimer(primer, seq):
    """
    Look for a primer sequence.

    @param primer: A C{str} primer sequence.
    @param seq: A BioPython C{Bio.Seq} sequence.

    @return: A C{list} of zero-based offsets into the sequence at which the
        primer can be found. If no instances are found, return an empty
        C{list}.
    """
    offsets = []
    seq = seq.upper()
    primer = primer.upper()
    primerLen = len(primer)
    discarded = 0
    offset = seq.find(primer)

    while offset > -1:
        offsets.append(discarded + offset)
        seq = seq[offset + primerLen :]
        discarded += offset + primerLen
        offset = seq.find(primer)

    return offsets


def findPrimerBidi(primer, seq):
    """
    Look for a primer in a sequence and its reverse complement.

    @param primer: A C{str} primer sequence.
    @param seq: A BioPython C{Bio.Seq} sequence.

    @return: A C{tuple} of two lists. The first contains (zero-based)
        ascending offsets into the sequence at which the primer can be found.
        The second is a similar list ascending offsets into the original
        sequence where the primer matches the reverse complemented of the
        sequence. If no instances are found, the corresponding list in the
        returned tuple must be empty.
    """
    # Note that we reverse complement the primer to find the reverse
    # matches. This is much simpler than reverse complementing the sequence
    # because it allows us to use findPrimer and to deal with overlapping
    # matches correctly.
    forward = findPrimer(primer, seq)
    reverse = findPrimer(reverse_complement(primer), seq)
    return forward, reverse


def findPrimerBidiLimits(primer, seq):
    """
    Report the extreme (inner) offsets of primer in a sequence and its
    reverse complement.

    @param primer: A C{str} primer sequence.
    @param seq: A BioPython C{Bio.Seq} sequence.

    @return: A C{tuple} of two C{int} offsets. The first is a (zero-based)
        offset into the sequence that is beyond the first instance (if any)
        of the primer. The second is an offset into the original sequence of
        the beginning of the last instance of the primer in the reverse
        complemented sequence.

        In other words, if you wanted to chop all instances of a primer
        out of a sequence from the start and the end (when reverse
        complemented) you'd call this function and do something like this:

          start, end = findPrimerBidiLimits(primer, seq)
          seq = seq[start:end]
    """
    forward, reverse = findPrimerBidi(primer, seq)
    if forward:
        start = forward[-1] + len(primer)
        end = len(seq)
        for offset in reverse:
            if offset >= start:
                end = offset
                break
    else:
        start = 0
        end = reverse[0] if reverse else len(seq)

    return start, end
