from unittest import TestCase

from dark import taxonomy


class BlastHits(object):
    """
    This is a placeholder to keep the test suite quiet now that we no
    longer have dark.blast.BlastHits. Taxonomy testing needs work.
    """


class FakeCursor(object):
    def __init__(self, results):
        self._results = results
        self._index = -1

    def execute(self, p):
        pass

    def fetchone(self):
        self._index += 1
        return self._results[self._index]

    def close(self):
        pass


class FakeDbConnection(object):
    def __init__(self, results):
        self._results = results
        self.open = True

    def cursor(self):
        return FakeCursor(self._results)

    def close(self):
        self.open = False


# TODO: put these taxonomy tests elsewhere after modernizing them.
class TestTaxonomy(object):
    """
    Test the helper functions in taxonomy.py
    """
    def testGetLineageInfo(self):
        """
        getLineageInfo should return the right taxIDs from
        the database
        """
        blastHits = BlastHits(None)
        blastHits.addHit('Smiley Cell Polyomavirus', {
            'taxID': 4,
        })
        db = FakeDbConnection([('species', 1), ('A')])
        result = taxonomy.getLineageInfo(blastHits, db=db)
        self.assertEqual({4: [{
                              'taxID': 4,
                              'parentTaxID': 1,
                              'rank': 'species',
                              'scientificName': 'A',
                              }]
                          }, result)

    def testGetLineageInfoDoesNotCloseOpenDatabase(self):
        """
        getLineageInfo should not close an open database it is passed.
        """
        blastHits = BlastHits(None)
        blastHits.addHit('Smiley Cell Polyomavirus', {
            'taxID': 4,
        })
        db = FakeDbConnection([('species', 1), ('A')])
        taxonomy.getLineageInfo(blastHits, db=db)
        self.assertEqual(True, db.open)

    def testTaxIDsPerTaxonomicRankEmptyInput(self):
        """
        taxIDsPerTaxonomicRank should return an empty dictionary
        when given an empty input.
        """
        result = taxonomy.taxIDsPerTaxonomicRank({}, 'species')
        self.assertEqual({}, result)

    def testTaxIDsPerTaxonomicRankAllTaxIDsPresent(self):
        """
        taxIDsPerTaxonomicRank should print the right taxID
        for the given taxonomic rank.
        """
        taxIDLookUpDict = {
            1: [{
                'taxID': 1,
                'parentTaxID': 4,
                'rank': 'species',
                'scientificName': 'mouse',
                }, {
                'taxID': 4,
                'parentTaxID': 7,
                'rank': 'genus',
                'scientificName': 'mouseian',
                }],
            2: [{
                'taxID': 2,
                'parentTaxID': 5,
                'rank': 'genus',
                'scientificName': 'dog',
                }],
            3: [{
                'taxID': 3,
                'parentTaxID': 6,
                'rank': 'family',
                'scientificName': 'cat',
                }]
        }
        blastHits = BlastHits(None)
        blastHits.addHit('Smiley Cell Polyomavirus', {
            'taxID': 1,
        })
        blastHits.addHit('Pink Sheep Virus', {
            'taxID': 2,
        })
        blastHits.addHit('Flying Elephant Making Virus', {
            'taxID': 3,
        })

        result = taxonomy.taxIDsPerTaxonomicRank(taxIDLookUpDict, 'species')
        self.assertEqual({'mouse': set([1])}, result)

    def testReadsPerTaxonomicRankEmptyInput(self):
        """
        readsPerTaxonomicRank should return an empty dictionary
        when given an empty input.
        """
        blastHits = BlastHits(None)
        result = taxonomy.readsPerTaxonomicRank({}, blastHits, 'species')
        self.assertEqual({}, result)

    def testReadsPerTaxonomicRank(self):
        """
        readsPerTaxonomicRank should print the right readNum
        for the given taxonomic rank.
        """
        taxIDLookUpDict = {
            1: [{
                'taxID': 1,
                'parentTaxID': 4,
                'rank': 'species',
                'scientificName': 'mouse',
                }, {
                'taxID': 4,
                'parentTaxID': 7,
                'rank': 'genus',
                'scientificName': 'mouseian',
                }],
            2: [{
                'taxID': 2,
                'parentTaxID': 5,
                'rank': 'genus',
                'scientificName': 'dog',
                }],
            3: [{
                'taxID': 3,
                'parentTaxID': 6,
                'rank': 'family',
                'scientificName': 'cat',
                }]
        }
        blastHits = BlastHits(None)
        blastHits.addHit('Smiley Cell Polyomavirus', {
            'taxID': 1,
            'plotInfo': {
                'items': [{
                    'readNum': 1234
                    }]
                }
        })
        blastHits.addHit('Pink Sheep Virus', {
            'taxID': 2,
            'plotInfo': {
                'items': [{
                    'readNum': 1235
                    }]
                }
        })
        blastHits.addHit('Flying Elephant Making Virus', {
            'taxID': 3,
            'plotInfo': {
                'items': [{
                    'readNum': 1236
                    }]
                }
        })

        result = taxonomy.readsPerTaxonomicRank(taxIDLookUpDict,
                                                blastHits, 'species')
        self.assertEqual({'mouse': set([1234])}, result)

    def testSubjectsPerTaxonomicRankEmptyInput(self):
        """
        subjectsPerTaxonomicRank should return an empty dictionary
        when given an empty input.
        """
        blastHits = BlastHits(None)
        result = taxonomy.subjectsPerTaxonomicRank({}, blastHits, 'species')
        self.assertEqual({}, result)

    def testSubjectsPerTaxonomicRank(self):
        """
        subjectsPerTaxonomicRank should print the right subject
        for the given taxonomic rank.
        """
        taxIDLookUpDict = {
            1: [{
                'taxID': 1,
                'parentTaxID': 4,
                'rank': 'species',
                'scientificName': 'mouse',
                }, {
                'taxID': 4,
                'parentTaxID': 7,
                'rank': 'genus',
                'scientificName': 'mouseian',
                }],
            2: [{
                'taxID': 2,
                'parentTaxID': 5,
                'rank': 'genus',
                'scientificName': 'dog',
                }],
            3: [{
                'taxID': 3,
                'parentTaxID': 6,
                'rank': 'family',
                'scientificName': 'cat',
                }]
        }
        blastHits = BlastHits(None)
        blastHits.addHit('Smiley Cell Polyomavirus', {
            'taxID': 1,
        })
        blastHits.addHit('Pink Sheep Virus', {
            'taxID': 2,
        })
        blastHits.addHit('Flying Elephant Making Virus', {
            'taxID': 3,
        })
        result = taxonomy.subjectsPerTaxonomicRank(taxIDLookUpDict,
                                                   blastHits, 'species')
        self.assertEqual({'mouse': set(['Smiley Cell Polyomavirus'])}, result)


class TestLineageFetcher(TestCase):
    """
    Test LineageFetcher class.
    """
    def testGetTaxonomy(self):
        pass
        # title = 'gi|5|gb|EU375804.1| Merkel cell polyomavirus'
        # lineageFetcher = taxonomy.LineageFetcher()
        #
        # db = FakeDbConnection([[15], [2], ['Merkel cell polyomavirus'],
        #                       [3], ['Polyomavirus'], [2],
        #                       ['dsDNA viruses'], [1], ['Vira']])
        # cursor = db.cursor()
        # lineage = getTaxonomy(title, cursor)
        # self.assertEqual(['Merkel cell polyomavirus', 'Polyomavirus',
        #                   'dsDNA viruses', 'Vira'], lineage)
