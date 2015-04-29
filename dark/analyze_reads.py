from Bio import SeqIO


def _longestPrefixOfTwoSeqs(a, b):
    length = min(len(a), len(b))
    result = 0
    while result < length:
        if a[result] == b[result]:
            result += 1
        else:
            break
    return result


def getPrefixAndSuffix(file_handle):
    read_list = list(SeqIO.parse(file_handle, 'fasta'))
    reversed_read_list = [read[::-1] for read in read_list]

    def longestCommonPrefix(read_list):
        sequences = read_list
        nSequences = len(sequences)
        if nSequences == 1:
            return len(sequences[0])
        elif nSequences == 0:
            return 0
        else:
            prefix = sequences[0]
            result = len(prefix)
            index = 1
            while index < nSequences:
                thisLen = _longestPrefixOfTwoSeqs(prefix, sequences[index])
                if thisLen == 0:
                    return 0
                elif thisLen < result:
                    prefix = prefix[:thisLen]
                    result = thisLen
                index += 1
            return result

    prefix = longestCommonPrefix(read_list)
    suffix = longestCommonPrefix(reversed_read_list)
    return prefix, suffix


def trimReads(prefix, suffix, file_handle):
    for record in SeqIO.parse(file_handle, 'fasta'):
        if suffix == 0:
            yield record[prefix:]
        else:
            yield record[prefix:-suffix]
