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


def getReads(filename):
    read_list = []
    for record in SeqIO.parse(open(filename, "rU"), "fasta"):
        read_list.append(record)
    print read_list

    reversed_read_list = map(lambda read: read[::-1], read_list)
    print reversed_read_list

    def longestCommonPrefix(read_list):
        sequences = read_list
        nSequences = len(sequences)
        if nSequences == 1:
            print "This file contains one read only!"
            return len(sequences[0])
        else:
            prefix = sequences[0]
            result = len(prefix)
            index = 1
            while index < nSequences:
                thisLen = _longestPrefixOfTwoSeqs(prefix, sequences[index])
                if thisLen == 0:
                    result = 0
                    return result
                elif thisLen < result:
                    prefix = prefix[:thisLen]
                    result = thisLen
                    return result
                index += 1

    result_prefix = longestCommonPrefix(read_list)
    result_suffix = longestCommonPrefix(reversed_read_list)
    print "result_prefix", result_prefix
    print "result_suffix", result_suffix

    def clipSequences(filename):
        print result_prefix
        print result_suffix
        for record in SeqIO.parse(open(filename, "rU"), "fasta"):
            if result_suffix == 0:
                print record[result_prefix:]
            else:
                print record[result_prefix:-result_suffix]
    clipSequences(filename)

