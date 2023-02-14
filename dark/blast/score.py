"""
Probably the most comprehensive online reference for the way BLAST
calculates scores is "BLAST Scoring Parameters" by E. Michael Gertz

You can get a copy at
    ftp://ftp.ncbi.nlm.nih.gov/blast/documents/developer/scoring.pdf

or
    https://notebooks.antigenic-cartography.org/terry/docs/papers/\
    blast-scoring-parameters.pdf

Section 5 has the calculation of effective search space.

For details on going between bit scores and e-values, see
The Statistics of Sequence Similarity Scores, at
http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
"""

from math import log

_LOG2 = log(2.0)


def bitScoreToEValue(bitScore, dbSize, dbSequenceCount, queryLength, lengthAdjustment):
    """
    Convert a bit score to an e-value.

    @param bitScore: The C{float} bit score to convert.
    @param dbSize: The C{int} total size of the database (i.e., the sum of
        the lengths of all sequences in the BLAST database).
    @param dbSequenceCount: The C{int} number of sequences in the database.
    @param queryLength: The C{int} length of the query.
    @param lengthAdjustment: The C{int} length adjustment (BLAST XML output
        calls this the Statistics_hsp-len).
    @return: A C{float} e-value.
    """
    effectiveDbSize = (dbSize - dbSequenceCount * lengthAdjustment) * (
        queryLength - lengthAdjustment
    )
    return effectiveDbSize * (2.0 ** (-1.0 * bitScore))


def eValueToBitScore(eValue, dbSize, dbSequenceCount, queryLength, lengthAdjustment):
    """
    Convert an e-value to a bit score.

    @param eValue: The C{float} e-value to convert.
    @param dbSize: The C{int} total size of the database (i.e., the sum of
        the lengths of all sequences in the BLAST database).
    @param dbSequenceCount: The C{int} number of sequences in the database.
    @param queryLength: The C{int} length of the query.
    @param lengthAdjustment: The C{int} length adjustment (BLAST XML output
        calls this the Statistics_hsp-len).
    @return: A C{float} bit score.
    """
    effectiveDbSize = (dbSize - dbSequenceCount * lengthAdjustment) * (
        queryLength - lengthAdjustment
    )
    return -1.0 * (log(eValue / effectiveDbSize) / _LOG2)
