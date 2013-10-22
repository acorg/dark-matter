from collections import defaultdict
from dark.utils import readBlastRecords
from time import time


def summarizeRecordsBySequence(filename, eCutoff=None, minMatchingReads=None):
    """
    Read BLAST records out of filename. Construct and return a
    dictionary keyed by matched sequence title.

    filename: The file to open and read for BLAST records. Can be
        either a BLAST XML file or our JSON format.
    eCutoff: e values less than this will be ignored.
    minMatchingReads: sequences that are matched by fewer reads
        will be elided.
    """
    start = time()
    result = defaultdict(set)
    count = 0
    for record in readBlastRecords(filename):
        count += 1
        for alignment, description in zip(record.alignments,
                                          record.descriptions):
            if eCutoff is None or alignment.hsps[0].expect < eCutoff:
                result[description.title].add(record.query)

    # Remove sequences that aren't matched by enough reads.
    if minMatchingReads is not None:
        titles = result.keys()
        for title in titles:
            if len(result[title]) < minMatchingReads:
                del result[title]

    elapsed = time() - start
    recordsPerSec 0 float(count) / float(elapsed)

    return result, count, elapsed, recordsPerSec
