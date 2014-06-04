from unittest import TestCase
from cStringIO import StringIO

from dark import fastamanipulations



class TestFastaSubset(TestCase):
    """
    Test the fastaSubset() function.
    """
    fasta = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
        atggctattgaactgtatct'

    def testNoReadIdsPresentTrue(self):
        readIds = []
        result = fastamanipulations.fastaSubset(StringIO(fasta), readIds, present=True)
        self.assertEqual(result, 1)

    def testNoReadIdsPresentFalse(self):
        readIds = []
        result = fastamanipulations.fastaSubset(StringIO(fasta), readIds, present=False)
        self.assertEqual(result, [])

    def testNoReadIdsInFasta(self):
        readIds = ['hola', 'there']
        result = fastamanipulations.fastaSubset(StringIO(fasta), readIds, present=True)
        self.assertEqual(result, 1)

    def testReadIdsPresentTrue(self):
        readIds = ['hey', 'you']
        result = fastamanipulations.fastaSubset(StringIO(fasta), readIds, present=True)
        self.assertEqual(result, '>hey\nagtcagtcagtc\n>you\nacctg')

    def testReadIdsPresentFalse(self):
        readIds = ['hey', 'you']
        result = fastamanipulations.fastaSubset(StringIO(fasta), readIds, present=True)
        self.assertEqual(result, 'how\natgggtc\n>are\natggctattgaactgtatct')


class TestFastaSubtract(TestCase):
    """
    Test the fastaSubtract() function.
    """
    fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
        atggctattgaactgtatct'
    fasta2 = '>all\nagtcagtcagtc\n>good\nacctg\n>over\natgggtc\n>here\n\
        atggctattgaactgtatct'
    fasta3 = '>hey\nagtcagtcagtc\n>you\nacctg\n>this\natgggtc\n>is\n\
        atggctattgaactgtatct'
    fasta4 = '>hey\nagtcagtcagtc\n>you\nacctg\n>nothing\natgggtc\n>there\n\
        atggctattgaactgtatct'
    fasta5 = '>this\nagtcagtcagtc\n>has\nacctg\n>many\natgggtc\n>rabbits\n\
        atggctattgaactgtatct'
    
    def testTwoFastaFiles(self):
        result = fastamanipulations.fastaSubtract([StringIO(file1), StringIO(file3)])
        self.assertEqual(result, ['hey', 'you'])

    def testThreeFastaFiles(self):
        result = fastamanipulations.fastaSubtract([StringIO(file1), StringIO(file3), StringIO(file4)])
        self.assertEqual(result, ['hey', 'you'])     

    def testTwoFastaFilesNoCommons(self):
        result = fastamanipulations.fastaSubtract([StringIO(file1), StringIO(file2)])
        self.assertEqual(result, []) 

    def testThreeFastaFilesNoCommons(self):
        result = fastamanipulations.fastaSubtract([StringIO(file1), StringIO(file2), StringIO(file5)])
        self.assertEqual(result, [])


class TestGetReadsIdsFromBlast(TestCase):
    """
    Test the getReadsIdsFromBlast() function.
    """
    params = {
            'application': 'BLASTN',
        }

        record = {
            'query': 'readEValue',
            'alignments': [
                {
                    'length': 37108,
                    'hsps': [
                        {
                            'bits': 20,
                            'sbjct_end': 15400,
                            'expect': 1e-100,
                            'sbjct': 'TACCC--CGGCCCGCG-CGGCCGGCTCTCCA',
                            'sbjct_start': 15362,
                            'query': 'TACCCTGCGGCCCGCTACGGCTGG-TCTCCA',
                            'frame': [1, 1],
                            'query_end': 68,
                            'query_start': 28
                        }
                    ],
                    'title': 'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                }
            ]
        }
        record = {
            'query': 'readBit',
            'alignments': [
                {
                    'length': 37108,
                    'hsps': [
                        {
                            'bits': 200,
                            'sbjct_end': 15400,
                            'expect': 1,
                            'sbjct': 'TACCC--CGGCCCGCG-CGGCCGGCTCTCCA',
                            'sbjct_start': 15362,
                            'query': 'TACCCTGCGGCCCGCTACGGCTGG-TCTCCA',
                            'frame': [1, 1],
                            'query_end': 68,
                            'query_start': 28
                        }
                    ],
                    'title': 'gi|887699|gb|DQ37780 Squirrelpox virus 1296/99',
                }
            ]
        }
        mockOpener = mockOpen(
            read_data=dumps(params) + '\n' + dumps(record) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            blastRecords = BlastRecords('file.json', None, None)
            self.assertEqual(1, len(list(blastRecords.records())))

    def testAreAllReadsBelowECutoffFound(self):
        with patch('__builtin__.open', mockOpener, create=True):
            readIds = fastamanipulations.getReadsIdsFromBlast('file.json', eCutoff=1e-2, bitCutoff=None)
            self.assertEqual(readids, ['readEValue'])

    def testAreAllReadsAboveBitCutoffFound(self):
        with patch('__builtin__.open', mockOpener, create=True):
            readIds = fastamanipulations.getReadsIdsFromBlast('file.json', eCutoff=None, bitCutoff=50)
            self.assertEqual(readids, ['readBit'])

