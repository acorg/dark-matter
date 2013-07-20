from Bio import SeqIO
from collections import defaultdict


def summarize_reads(filename, returnSequences=False):
    """
    open a .fasta file, prints number of of reads, average length of read,
    total number of bases, longest, shortest and median read, total number
    and average of individual base (A, T, G, C, N)
    """
    base_counts = defaultdict(int)
    read_number = 0
    total_length = 0
    length_list = []
    sequenceList = []

    for record in SeqIO.parse(open(filename, "rU"), "fasta"):
        total_length += len(record)
        read_number += 1
        length_list.append(len(record))
        sequenceList.append(record)
        for base in record:
            base_counts[base] += 1

    def median(length_list):
        my_list = sorted(length_list)
        list_count = len(length_list)
        if list_count % 2 != 0:
            index_middle_number = (list_count - 1) / 2
            return my_list[index_middle_number]
        else:
            first_number = my_list[list_count/2]
            second_number = my_list[(list_count/2)-1]
            return (first_number+second_number)/2.0

    result = {
        "read_number": read_number,
        "total_length": total_length,
        "average_length": total_length / read_number,
        "max_length": max(length_list),
        "min_length": min(length_list),
        "median_length": median(length_list),
        "base_counts": base_counts
    }

    if returnSequences:
        result["sequences"] = sequenceList

    return result
