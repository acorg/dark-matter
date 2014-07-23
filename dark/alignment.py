class Alignment(object):
    """
    Hold information about a read alignment.

    @param hitLength: The C{int} length of the sequence the read hit against.
    @param hitTitle: The C{str} title of the sequence the read hit against.
    """

    def __init__(self, hitLength, hitTitle):
        self.hitLength = hitLength
        self.hitTitle = hitTitle
        self.hsps = []

    def addHsp(self, hsp):
        """
        Add an HSP to the list of HSPs for this alignment.

        @param hsp: A L{dark.hsp} (or subclass) instance.
        """
        self.hsps.append(hsp)
