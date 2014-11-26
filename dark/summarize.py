from Bio import SeqIO
from collections import defaultdict, Counter


def summarizeReads(file_handle, file_type):
    """
    open a fasta or fastq file, prints number of of reads,
    average length of read, total number of bases, longest,
    shortest and median read, total number and average of
    individual base (A, T, G, C, N).
    """
    base_counts = defaultdict(int)
    read_number = 0
    total_length = 0
    length_list = []

    records = SeqIO.parse(file_handle, file_type)

    for record in records:
        total_length += len(record)
        read_number += 1
        length_list.append(len(record))
        for base in record:
            base_counts[base] += 1

    def median(length_list):
        my_list = sorted(length_list)
        list_count = len(length_list)
        if list_count % 2 != 0:
            index_middle_number = (list_count - 1) / 2
            return my_list[index_middle_number]
        else:
            first_number = my_list[list_count / 2]
            second_number = my_list[(list_count / 2) - 1]
            return (first_number + second_number) / 2.0

    result = {
        "read_number": read_number,
        "total_length": total_length,
        "average_length": total_length / read_number if read_number > 0 else 0,
        "max_length": max(length_list) if length_list else 0,
        "min_length": min(length_list) if length_list else 0,
        "median_length": median(length_list) if length_list else 0,
        "base_counts": base_counts
    }

    return result


def summarizePosition(records, position):
    """
    For a fasta file with sequences, summarize what is happening at a specific
    position.

    @param fileName: a C{str} fileName of a fasta file.
    @param position: a position in the sequence that should be looked at. This
        assumes that the position is given in 1-based offsets (NOT how python
        usually does it!).
    """
    countAtPosition = Counter()

    bases = []
    for record in records:
        sequence = record.sequence
        try:
            bases.append(sequence[position - 1])
        except IndexError:
            continue

    for base in bases:
        countAtPosition[base] += 1

    result = {
        'includedSequences': len(bases),
        'allSequences': len(records),
        'countAtPosition': countAtPosition
    }

    return result
