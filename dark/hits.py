def bestAlignment(hit):
    """
    Find the best alignment for a read hit. This is the one whose
    first HSP has the best score.

    Note that the comparison of HSP score values is taken care of by
    the HSP class. This works whether higher or lower scores are
    considered better.

    @param hit: a C{Hit} instance.
    @return: The alignment with the best first HSP score.
    """
    return max(hit.alignments,
               key=lambda alignment: alignment.hsps[0])


class Hit(object):
    """
    Holds information about a single read hit.

    @param read: A C{Read} instance.
    @param alignments: A C{list} of alignments. Each alignment is a C{dict} of
        the following form:

        {
          "hsps": [
            {
              "score": 293.432,
              "frame": [ 1, 1 ],  # Optional, used if nucleotides are involved.
              "readMatchedSequence": "CGCC",
              "readEnd": 492,
              "readStart": 31,
              "readEndInHit": 2926446,
              "readStartInHit": 2925990,
              "hitMatchedSequence": "CGCC",
              "hitEnd": 2926446,
              "hitStart": 2925990
            },
            # Optionally, additional hsp dicts.
          ],
          "hitLength": 5371077,
          "hitTitle": "gi|2577609|dbj|AP01960.1| Escherichia coli DNA"
        }

    """
    def __init__(self, read, alignments):
        self.read = read
        self.alignments = alignments


class Hits(object):
    """
    Maintain a collection of read hits.

    @param reads: A L{Reads} instance, containing the reads (sequences)
        given to the application to create these hits.
    @param application: The C{str} name of the application that generated
        these read hits.
    @param params: A C{dict} of the parameters that were given to the
        application to create these hits.
    """

    def __init__(self, reads, application, params):
        self.reads = reads
        self.application = application
        self.params = params

    def compareScores(self, score1, score2):
        """
        Compare two HSP scores.

        @param score1: A numeric score.
        @param score2: A numeric score.
        @return: An C{int} that is negative, zero, or positive according to
            whether score1 is less than, equal to, or greater than score2.
        """
        return cmp(score1, score2)

    def filter(self, limit=None, oneAlignmentPerRead=False, maxHspsPerHit=None,
               scoreCutoff=None):
        """
        @param limit: An C{int} limit on the number of records to read.
        @param oneAlignmentPerRead: if C{True}, only keep the best
            alignment for each read.
        @param maxHspsPerHit: The maximum number of HSPs to keep for each
            alignment for each read.
        @param scoreCutoff: A C{float} score. Hits with scores that are not
            better than this score will be ignored.
        @return: A generator that yields L{Hit} instances.
        """

        count = 0
        for hit in self:
            if limit is not None and count == limit:
                return

            # Throw out any unwanted HSPs due to maxHspsPerHit.
            for alignment in hit.alignments:
                hsps = alignment.hsps
                if maxHspsPerHit is not None and len(hsps) > maxHspsPerHit:
                    alignment.hsps = hsps[:maxHspsPerHit]

            if oneAlignmentPerRead and hit.alignments:
                hit.alignments = [bestAlignment(hit)]

            if scoreCutoff is not None:
                wantedAlignments = []
                for alignment in hit.alignments:
                    hsps = alignment.hsps
                    wantedHsps = []
                    for hsp in hsps:
                        if hsp.betterThan(scoreCutoff):
                            wantedHsps.append(hsp)
                    if wantedHsps:
                        alignment.hsps = wantedHsps
                        wantedAlignments.append(alignment)

                if wantedAlignments:
                    hit.alignments = wantedAlignments
                else:
                    continue

            yield hit
            count += 1
