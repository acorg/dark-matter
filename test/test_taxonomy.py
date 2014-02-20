from unittest import TestCase

from dark.blast import BlastHits
from dark import taxonomy


class FakeCursor(object):
    result = [['species', 3],
               'Smileyus viriae',
              ['genus', 4], 'Lucky viriae']
    def __init__(self, result):
        self._result = result
        self._index = -1

    def execute(self):
        self._index += 1
        return self._result[self._index]

    def close():
        pass


class FakeDbConnection(object):
    def cursor(self):
        db = FakeCursor(object)
        return db

    def close(self):
        pass


class TestTaxonomy(TestCase):
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
            'taxID': 1,
        })
        blastHits.addHit('Pink Sheep Virus', {
            'taxID': 2,
        })
        db = FakeDbConnection()
        result = taxonomy.getLineageInfo(blastHits, db=db)
        self.assertEqual({1: [{
                              'taxID': 1,
                              'parentTaxID': 3,
                              'rank': 'species',
                              'scientificName': 'Smileyus viriae',
                              }],
                          2: [{'taxID': 2,
                               'parentTaxID': 4,
                               'rank': 'genus',
                               'scientificName': 'Lucky viriae'
                               }]
                          }, result)

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
