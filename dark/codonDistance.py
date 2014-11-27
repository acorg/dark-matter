from collections import defaultdict


def findDistance(co1, co2):
    """
    Find the distance between two codons.
    """
    inequal = 0
    count = 0
    while count < 3:
        if co1[count] != co2[count]:
            inequal += 1
        count += 1
    return inequal


def codonInformation(codons1, codons2):
    """
    Take two C{list} of codons, and returns information about a and the min and
    max distance between them.

    @param codon1: a C{list} of codons.
    @param codon2: a C{list} of codons.

    @return: a dict with how often a distance occurs between which codons.
    """
    result = defaultdict(list)
    for c1 in codons1:
        for c2 in codons2:
            distance = findDistance(c1, c2)
            result[distance].append([c1, c2])

    return result
