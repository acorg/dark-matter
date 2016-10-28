from time import time
from collections import defaultdict

# TODO: Some of these imported functions are no longer in utils.py!
from dark.utils import findHits, getSequence, summarizeHits, printHSP, report


def consensusSequence(recordFilename, hitId, fastaFilename, eCutoff=None,
                      db='nt', actualSequenceId=None):
    """
    Build a consensus sequence against a target sequence.

    recordFilename: the BLAST XML output file.
    hitId: the str sequence id to examine the BLAST output for hits against.
    fastaFilename: The name of the FASTA file containing query sequences,
        or a list of fasta sequences from SeqIO.parse
    eCutoff: converted e values less than this will be ignored.
    db: the BLAST db to use to look up the target and, if given, actual
        sequence.
    actualSequenceId: the str id of the actual sequence (if known).
    """

    print('TODO: This function is not finished yet.')
    return

    start = time()
    if isinstance(recordFilename, str):
        # TODO: REMOVE THE LIMIT IN THE NEXT LINE!
        allhits = findHits(recordFilename, set([hitId]), limit=100)
    else:
        allhits = recordFilename
    sequence = getSequence(hitId, db)
    fasta, summary = summarizeHits(allhits, fastaFilename, eCutoff=eCutoff)
    minX, maxX = summary['minX'], summary['maxX']
    if actualSequenceId:
        # UNUSED.
        # actualSequence = getSequence(actualSequenceId, db)
        pass
    print(summary['hitCount'])
    print('seq len =', len(sequence))
    fasta = summary['fasta']
    # The length of the consensus depends on where the query sequences fell
    # when aligned with the target. The consensus could extend the target
    # at both ends.
    consensusLen = maxX - minX
    consensus = [None, ] * consensusLen
    for item in summary['items']:
        print('NEW HSP')
        printHSP(item['origHsp'])  # TODO: REMOVE ME
        hsp = item['hsp']
        print('HIT query-start=%d query-stop=%d subj-start=%d subj-stop=%d' % (
            hsp['queryStart'], hsp['queryEnd'], hsp['subjectStart'],
            hsp['subjectEnd']))
        # print '   match: %s%s' % ('.' * hsp['subjectStart'], '-' *
        # (hsp['subjectEnd'] - hsp['subjectStart']))
        if item['frame']['subject'] > 0:
            query = fasta[item['sequenceId']].seq
        else:
            query = fasta[item['sequenceId']].reverse_complement().seq
        print('       target:',
              sequence[hsp['queryStart']:hsp['queryEnd']].seq)
        print('        query:', query)
        match = []
        for index in range(hsp['subjectStart'], hsp['subjectEnd']):
            queryIndex = index - hsp['queryStart']
            match.append('.' if query[queryIndex] == sequence[index] else '*')
        print('        match: %s%s' % (
            ' ' * (hsp['subjectStart'] - hsp['queryStart']), ''.join(match)))
        print('        score:', item['convertedE'])
        print('    match len:', hsp['subjectEnd'] - hsp['subjectStart'])
        print('subject frame:', item['frame']['subject'])
        for queryIndex, sequenceIndex in enumerate(
                range(hsp['queryStart'], hsp['queryEnd'])):
            consensusIndex = sequenceIndex + minX
            locus = consensus[consensusIndex]
            if locus is None:
                consensus[consensusIndex] = locus = defaultdict(int)
            locus[query[queryIndex]] += 1

    # Print the consensus before the target, if any.
    for index in range(minX, 0):
        consensusIndex = index - minX
        if consensus[consensusIndex]:
            print('%d: %r' % (index, consensus[consensusIndex]))
    # Print the consensus as it overlaps with the target, if any.
    for index in range(0, len(sequence)):
        consensusIndex = index - minX
        try:
            if consensus[consensusIndex]:
                print('%d: %r (%s)' % (
                    index, consensus[index], sequence[index]))
        except KeyError:
            # There's nothing left in the consensus, so we're done.
            break
    for index in range(len(sequence), maxX):
        consensusIndex = index - minX
        if consensus[consensusIndex]:
            print('%d: %r' % (index, consensus[consensusIndex]))
    stop = time()
    report('Consensus sequence generated in %.3f mins.' %
           ((stop - start) / 60.0))
    return summary, consensus
