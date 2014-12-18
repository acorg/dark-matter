from collections import defaultdict


def findDistance(co1, co2):
    """
    Find the distance between two codons.

    @param co1: A C{str} of length three.
    @param co2: A C{str} of length three.
    @return: An C{int} distance in [0, 3].
    """
    return (co1[0] != co2[0]) + (co1[1] != co2[1]) + (co1[2] != co2[2])


def codonInformation(codons1, codons2):
    """
    Take two C{list} of codons, and returns information about the min and
    max distance between them.

    @param codon1: a C{list} of codons.
    @param codon2: a C{list} of codons.

    @return: a dict whose keys are C{int} distances and whose values are lists
        of pairs of codons.
    """
    result = defaultdict(list)
    for c1 in codons1:
        for c2 in codons2:
            distance = findDistance(c1, c2)
            result[distance].append([c1, c2])

    return result
