from collections import defaultdict


def titleCounts(readsAlignments):
    """
    Count the number of times each title in a readsAlignments instance is
    matched. This is useful for rapidly discovering what titles were matched
    and with what frequency.

    @param readsAlignments: A L{dark.alignments.ReadsAlignments} instance.
    @return: A C{dict} whose keys are titles and whose values are the integer
        counts of the number of reads that matched that title.
    """

    titles = defaultdict(int)
    for readAlignments in readsAlignments:
        for alignment in readAlignments.alignments:
            titles[alignment.subjectTitle] += 1
    return titles


class TitleAlignment(object):
    """
    Hold information about a read's HSPs for a title alignment.

    @param read: The C{Read} that aligned.
    @param hsps: A C{list} of L{dark.hsp.HSP} (or subclass) instance.
    """

    def __init__(self, read, hsps):
        self.read = read
        self.hsps = hsps


class TitleAlignments(object):
    """
    Holds information about a set of alignments against a sequence.

    @param subjectTitle: The C{str} title of the sequence the read hit against.
    @param subjectLength: The C{int} length of the sequence the read hit
        against.
    """

    def __init__(self, subjectTitle, subjectLength):
        # TODO: Do we need the title in here?
        self.subjectTitle = subjectTitle
        self.subjectLength = subjectLength
        self.alignments = []

    def addAlignment(self, alignment):
        """
        Add an alignment to the list of alignments that matched this title.

        @param alignment: A L{TitleAlignment} instance.
        """
        self.alignments.append(alignment)


class TitlesAlignments(dict):
    """
    Holds (as a dictionary) a set of titles, each with its alignments.

    @param readsAlignments: A L{dark.alignments.ReadsAlignments} instance.
    """

    def __init__(self, readsAlignments):
        dict.__init__(self)
        for readAlignments in readsAlignments:
            for alignment in readAlignments.alignments:
                title = alignment.subjectTitle
                try:
                    titleAlignments = self[title]
                except KeyError:
                    titleAlignments = self[title] = TitleAlignments(
                        title, alignment.subjectLength)
                titleAlignments.addAlignment(
                    TitleAlignment(readAlignments.read, alignment.hsps))
