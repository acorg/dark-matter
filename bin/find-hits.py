#!/usr/bin/env python

from dark.utils import fastaToList, findHits, getSequence, printHSP
import sys

if len(sys.argv) != 4:
    print >>sys.stderr, 'Usage: %s record-file.xml hit-id, fastaFile' % sys.argv[0]
    sys.exit(1)
else:
    recordFile, hitId, fastaFile = sys.argv[1:]
    target = str(getSequence(hitId).seq)
    fasta = [str(s.seq) for s in fastaToList(fastaFile)]
    for sequenceNum, hsps in findHits(recordFile, hitId):
        print 'Sequence', sequenceNum
        print fasta[sequenceNum]
        for i, hsp in enumerate(hsps, start=1):
            print '  HSP %d' % i
            printHSP(hsp, '    ')
            # print '    query ', fasta[sequenceNum][hsp.query_start - 1:hsp.query_end]
            # print '          ', hsp.match
            # print '    target', target[hsp.sbjct_start - 1:hsp.sbjct_end]
