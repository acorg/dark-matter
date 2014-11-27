import sys

from dark.codonTable import CODONS


def findDistance(codon1, codon2):
    """
    Find the distance between two codons.
    """
    inequal = 0
    count = 0
    while count < 3:
        if codon1[count] != codon2[count]:
            inequal += 1
        count += 1
    return inequal


def codonInformation(aa1, aa2):
    """
    Take two amino acids, in three letter code, and return information about
    which codons code for that aa and the min and max distance between them.
    """
    result = {}
    try:
        codonsAA1 = CODONS[aa1.lower()]
        codonsAA2 = CODONS[aa2.lower()]
    except KeyError:
        print >>sys.stderr, ('You must give you amino acids in three-letter '
                             'code.')
        sys.exit(1)

    result = {'codons': {aa1: codonsAA1,
                         aa2: codonsAA2},
              'distances': {0: [],
                            1: [],
                            2: [],
                            3: []
                            }
              }

    for c1 in codonsAA1:
        for c2 in codonsAA2:
            distance = findDistance(c1, c2)
            result['distances'][distance].append([c1, c2])

    return result
