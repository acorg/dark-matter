from typing import Union, Optional
from io import TextIOWrapper
from collections import defaultdict
from Bio import SeqIO

from dark.reads import Read
from dark.utils import median


def summarizeReads(
    file_handle: Union[str, TextIOWrapper], file_type: str
) -> dict[str, Union[int, float, defaultdict[str, int]]]:
    """
    open a fasta or fastq file, prints number of of reads,
    average length of read, total number of bases, longest,
    shortest and median read, total number and average of
    individual base (A, T, G, C, N).
    """
    base_counts: defaultdict[str, int] = defaultdict(int)
    read_number = 0
    total_length = 0
    length_list: list[int] = []

    records = SeqIO.parse(file_handle, file_type)

    for record in records:
        total_length += len(record)
        read_number += 1
        length_list.append(len(record))
        for base in record:
            base_counts[base] += 1

    result: dict[str, Union[int, float, defaultdict[str, int]]] = {
        "read_number": read_number,
        "total_length": total_length,
        "average_length": total_length / read_number if read_number > 0 else 0,
        "max_length": max(length_list) if length_list else 0,
        "min_length": min(length_list) if length_list else 0,
        "median_length": median(length_list) if length_list else 0,
        "base_counts": base_counts,
    }

    return result


def sequenceCategoryLengths(
    read: Read,
    categories: dict[str, str],
    defaultCategory: Optional[str] = None,
    suppressedCategory: str = "...",
    minLength: int = 1,
):
    """
    Summarize the nucleotides or AAs found in a read by assigning each to a
    category and reporting the lengths of the contiguous category classes
    found along the sequence.

    @param read: A C{Read} instance or one of its subclasses.
    @param categories: A C{dict} mapping nucleotides or AAs to category.
    @param defaultCategory: The category to use if a sequence base is not
        in C{categories}.
    @param suppressedCategory: The category to use to indicate suppressed
        sequence regions (i.e., made up of stretches of bases that are less
        than C{minLength} in length).
    @param minLength: stretches of the read that are less than this C{int}
        length will be summed and reported as being in the
        C{suppressedCategory} category.
    @raise ValueError: If minLength is less than one.
    @return: A C{list} of 2-C{tuples}. Each tuple contains a (category, count)
        where category is either a C{str} or C{None} and count is an C{int}.
    """
    result: list[tuple[Optional[str], int]] = []
    append = result.append
    get = categories.get
    first = True
    currentCategory = None
    currentCount = 0
    suppressing = False
    suppressedCount = 0

    if minLength < 1:
        raise ValueError("minLength must be at least 1")

    for base in read.sequence:
        thisCategory = get(base, defaultCategory)
        if first:
            first = False
            currentCategory = thisCategory
            currentCount += 1
        else:
            if thisCategory == currentCategory:
                # This base is still in the same category as the last base.
                # Keep counting.
                currentCount += 1
            else:
                # This is a new category.
                if currentCount < minLength:
                    # The category region that was just seen will not be
                    # emitted.
                    if suppressing:
                        # Already suppressing. Suppress the just-seen region too.
                        suppressedCount += currentCount
                    else:
                        # Start suppressing.
                        suppressedCount = currentCount
                        suppressing = True
                else:
                    if suppressing:
                        append((suppressedCategory, suppressedCount))
                        suppressedCount = 0
                        suppressing = False
                    append((currentCategory, currentCount))
                currentCategory = thisCategory
                currentCount = 1

    if suppressing:
        append((suppressedCategory, suppressedCount + currentCount))
    elif currentCount >= minLength:
        append((currentCategory, currentCount))
    elif currentCount:
        append((suppressedCategory, currentCount))

    return result
