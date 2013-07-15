from Bio import SeqIO
from hashlib import md5


def fastaToList(fastaFilename):
    return list(SeqIO.parse(fastaFilename, 'fasta'))


def dedupFasta(sequences):
    """
    sequences: an iterator producing Bio.Seq sequences.

    return: a generator of sequences with no duplicates.
    """
    seen = set()
    for s in sequences:
        thisSeq = str(s.seq)
        if thisSeq not in seen:
            seen.add(thisSeq)
            yield s


def dePrefixAndSuffixFasta(sequences):
    """
    sequences: an iterator producing Bio.Seq sequences.

    return: a generator of sequences with no duplicates and no fully contained
        subsequences.
    """
    sequences = sorted(sequences, key=lambda s: len(s.seq), reverse=True)
    seen = set()
    for s in sequences:
        thisSeq = str(s.seq)
        thisHash = md5(thisSeq).digest()
        if thisHash not in seen:
            # Add prefixes.
            newHash = md5()
            for nucl in thisSeq:
                newHash.update(nucl)
                seen.add(newHash.digest())
            # Add suffixes.
            for start in xrange(len(thisSeq) - 1):
                seen.add(md5(thisSeq[start + 1:]).digest())
            yield s
    # print 'seen contains %d hashes. total bytes = %d' % (len(seen),
    # len(seen) * len(md5('x').digest()))
