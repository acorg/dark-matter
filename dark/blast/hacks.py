from __future__ import print_function


def printBlastRecordForDerek(record):
    """
    This is a hacked version of printBlastRecord that I wrote in
    Addenbrookes hospital so Derek could try hacking on the reads
    that had low e values.  It's called by bin/print-blast-xml-for-derek.py
    """
    if len(record.alignments) == 0:
        return 0
    if record.descriptions[0].e > 1e-120:
        return 0
    print(record.query, end=' ')
    for i, alignment in enumerate(record.alignments):
        print(record.descriptions[i].e, end=' ')
        for hspIndex, hsp in enumerate(alignment.hsps, start=1):
            if hsp.sbjct_start < hsp.sbjct_end:
                sense = 1
                if hsp.sbjct_start < 2200:
                    side = 'left'
                else:
                    side = 'right'
            else:
                sense = -1
                if hsp.sbjct_end < 2200:
                    side = 'left'
                else:
                    side = 'right'
            print(sense, side)
            break
    return 1
