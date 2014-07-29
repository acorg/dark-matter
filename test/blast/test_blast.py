from dark.blast.blast import BlastHits


class TestTaxonomyFiltering:
    def testAllTaxonomy(self):
        blastHits = BlastHits(None)
        blastHits.addHit('gi|293595919|gb|HM011539.1|', {
            'eMin': 0.01,
            'taxonomy': ['Merkel cell polyomavirus',
                         'unclassified Polyomavirus',
                         'Polyomavirus', 'Polyomaviridae',
                         'dsDNA viruses, no RNA stage', 'Vira']
        })
        result = blastHits.filterHits(taxonomy='all')
        self.assertEqual(result.titles, {'gi|293595919|gb|HM011539.1|': {
            'eMin': 0.01,
            'taxonomy': ['Merkel cell polyomavirus',
                         'unclassified Polyomavirus',
                         'Polyomavirus', 'Polyomaviridae',
                         'dsDNA viruses, no RNA stage',
                         'Vira']}})

    def testGivenTaxonomyNotPresent(self):
        blastHits = BlastHits(None)
        blastHits.addHit('gi|293595919|gb|HM011539.1|', {
            'eMin': 0.01,
            'taxonomy': ['Merkel cell polyomavirus',
                         'unclassified Polyomavirus',
                         'Polyomavirus', 'Polyomaviridae',
                         'dsDNA viruses, no RNA stage',
                         'Vira']
        })
        result = blastHits.filterHits(taxonomy='fiction')
        self.assertEqual(result.titles, {})

    def testGivenTaxonomyPresent(self):
        blastHits = BlastHits(None)
        blastHits.addHit('gi|293595919|gb|HM011539.1|', {
            'eMin': 0.01,
            'taxonomy': ['Merkel cell polyomavirus',
                         'unclassified Polyomavirus',
                         'Polyomavirus', 'Polyomaviridae',
                         'dsDNA viruses, no RNA stage',
                         'Vira']
        })
        result = blastHits.filterHits(taxonomy='Vira')
        self.assertEqual(result.titles, {'gi|293595919|gb|HM011539.1|': {
                                         'eMin': 0.01,
                                         'taxonomy': ['Merkel cell '
                                                      'polyomavirus',
                                                      'unclassified '
                                                      'Polyomavirus',
                                                      'Polyomavirus',
                                                      'Polyomaviridae',
                                                      'dsDNA viruses, '
                                                      'no RNA stage',
                                                      'Vira']}})
