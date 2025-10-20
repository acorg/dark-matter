from typing import Iterator, Sequence

from dark import File
from dark.fasta import FastaReads
from dark.reads import Read


def _longestPrefixOfTwoSeqs(a: str, b: str) -> int:
    length = min(len(a), len(b))
    result = 0
    while result < length:
        if a[result] == b[result]:
            result += 1
        else:
            break
    return result


def getPrefixAndSuffix(file_handle: File) -> tuple[int, int]:
    read_list = list(FastaReads(file_handle))
    reversed_read_list = [read[::-1] for read in read_list]

    def longestCommonPrefix(read_list: Sequence) -> int:
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
                thisLen = _longestPrefixOfTwoSeqs(
                    prefix.sequence, sequences[index].sequence
                )
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


def trimReads(prefix: int, suffix: int, file_handle: File) -> Iterator[Read]:
    for record in FastaReads(file_handle):
        if suffix == 0:
            yield record[prefix:]
        else:
            yield record[prefix:-suffix]
