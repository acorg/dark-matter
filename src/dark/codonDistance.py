from collections import defaultdict

TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def findDistance(co1, co2):
    """
    Find the distance between two codons.

    @param co1: A C{str} of length three.
    @param co2: A C{str} of length three.
    @return: An C{int} distance in [0, 3].
    """
    return (co1[0] != co2[0]) + (co1[1] != co2[1]) + (co1[2] != co2[2])


def codonInformation(codons1, codons2, countTransitions=False):
    """
    Take two C{list} of codons, and returns information about the min and
    max distance between them.

    @param codon1: a C{list} of codons.
    @param codon2: a C{list} of codons.
    @param countTransitions: A C{bool}. If C{True}, the tuples in the C{list}
        values of the returned C{dict} will also contain a count of the number
        of transitions (as opposed to transversions) between the pair of codons.

    @return: a dict whose keys are C{int} distances and whose values are tuples
        of pairs of codons (and the number of transitions involved, if
        C{countTransitions} is C{True})
    """
    result = defaultdict(list)
    for c1 in codons1:
        for c2 in codons2:
            distance = findDistance(c1, c2)
            if countTransitions:
                transitions = sum((a, b) in TRANSITIONS for a, b in zip(c1, c2))
                result[distance].append((c1, c2, transitions))
            else:
                result[distance].append((c1, c2))

    return result
