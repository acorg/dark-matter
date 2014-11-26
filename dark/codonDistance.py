from dark.codonTable import CODONS


def codonInformation(aa1, aa2):
    """
    Take two amino acids, in three letter code, and return information about
    which codons code for that aa and the min and max distance between them.
    """
    result = {}

    codonsAA1 = CODONS[aa1]
    codonsAA2 = CODONS[aa2]

    result['codons'] = {aa1: codonsAA1,
                        aa2: codonsAA2}

    def findDistance(codon1, codon2):
        """
        Find the distance between two codons.
        """
        print codon1, codon2
        inequal = 0
        count = 0
        while count < 3:
            if codon1[count] != codon2[count]:
                inequal += 1
            count += 1
        return inequal

    result['minCodons'] = []
    result['maxCodons'] = []
    result['minDistance'] = 10
    result['maxDistance'] = 0
    for c1 in codonsAA1:
        for c2 in codonsAA2:
            distance = findDistance(c1, c2)
            if distance < result['minDistance']:
                result['minDistance'] = distance
                result['minCodons'] = [c1, c2]
            elif distance > result['maxDistance']:
                result['maxDistance'] = distance
                result['maxCodons'] = [c1, c2]
            elif distance == result['minCodons']:
                result['minDistance'].extend([c1, c2])
            elif distance == result['maxDistance']:
                result['maxCodons'].extend([c1, c2])
    print result

    return result
