from dark import mysql


class BlastHits(object):
    """
    Maintain information about a set of sequences hit by a BLAST run.

    @param records: A L{BlastRecords} instance.
    @param readSetFilter: A L{dark.filter.ReadSetFilter} instance or C{None}.
        If not C{None}, this is the read set filter that was used to filter
        the hits for this read set.
    """

    def __init__(self, records, readSetFilter=None):
        self.records = records
        self.readSetFilter = readSetFilter
        self.titles = {}
        self.fasta = None  # Computed (once) in summarizeHits.
        self.plotParams = None  # Set in computePlotInfo.

    def _convertBitScoresToRanks(self, plotInfo):
        """
        Change the bit scores for the reads that hit each sequence to be their
        ranks.

        @param plotInfo: A C{dict} of plot information, as returned by
            C{plotInfoDict} in C{computePlotInfo}.
        """
        items = plotInfo['items']
        items.sort(key=lambda item: item['bitScore'])
        for rank, item in enumerate(items, start=1):
            item['bitScore'] = rank
        plotInfo['bitScoreMax'] = len(items)
        plotInfo['bitScoreMin'] = 1

    def getTaxIDFromMySql(self, db=None):
        """
        For each title in C{self.titles}, read the corresponding taxId from
        the gi_taxid_nucl table in a MySQL database called ncbi_taxonomy.
        """
        if db is None:
            openedDb = True
            db = mysql.getDatabaseConnection()
        else:
            openedDb = False
        cursor = db.cursor()

        # for each title (=gi number) get the taxId from the database
        # and add it to self.titles.
        for title in self.titles:
            giNr = int(title.split('|')[1])
            question = 'SELECT taxID from gi_taxid_nucl where gi = %d' % giNr
            cursor.execute(question)
            try:
                result = cursor.fetchone()[0]
            except TypeError:
                result = 'No taxID found'
            self.titles[title]['taxID'] = result

        cursor.close()
        if openedDb:
            db.close()
