from unittest import TestCase
from cStringIO import StringIO
from mocking import mockOpen
from json import dumps
from mock import patch
from Bio import SeqIO, SeqRecord, Seq
from dark import fastamanipulations


class TestFastaSubset(TestCase):
    """
    Test the fastaSubset() function.
    """

    def testNoReadIdsPresentTrue(self):
        fasta = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        readIds = []
        result = fastamanipulations.fastaSubset(StringIO(fasta),
                                                readIds, present=True)
        self.assertEqual(result, None)

    def testNoReadIdsPresentFalse(self):
        fasta = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        readIds = []
        result = fastamanipulations.fastaSubset(StringIO(fasta),
                                                readIds, present=False)
        self.assertEqual(result, '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\natggctattgaactgtatct')

    def testNoReadIdsInFasta(self):
        fasta = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        readIds = ['hola', 'there']
        result = fastamanipulations.fastaSubset(StringIO(fasta),
                                                readIds, present=True)
        self.assertEqual(result, None)

    def testReadIdsPresentTrue(self):
        fasta = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fakeResult = SeqIO.read(StringIO('>are\natggctattgaactgtatct'), 'fasta')
        readIds = ['hey', 'you', 'how']
        result = fastamanipulations.fastaSubset(StringIO(fasta),
                                                readIds, present=True)
        self.assertEqual(result, [fakeResult])

    def testReadIdsPresentFalse(self):
        fasta = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fakeResult = SeqIO.read(StringIO('>hey\nagtcagtcagtc'), 'fasta')
        readIds = ['hey']
        result = fastamanipulations.fastaSubset(StringIO(fasta),
                                                readIds, present=True)
        self.assertEqual(result, [fakeResult])


class TestFastaSubtract(TestCase):
    """
    Test the fastaSubtract() function.
    """

    def testTwoFastaFiles(self):
        fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta3 = '>hey\nagtcagtcagtc\n>you\nacctg\n>this\natgggtc\n>is\n\
            atggctattgaactgtatct'

        result = fastamanipulations.fastaSubtract([StringIO(fasta1),
                                                  StringIO(fasta3)])
        self.assertEqual(result, [])

    def testThreeFastaFiles(self):
        fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta3 = '>hey\nagtcagtcagtc\n>you\nacctg\n>this\natgggtc\n>is\n\
            atggctattgaactgtatct'
        fasta4 = '>hey\nagtcagtcagtc\n>you\nacctg\n>nothing\natgggtc\n>there\n\
            atggctattgaactgtatct'

        result = fastamanipulations.fastaSubtract([StringIO(fasta1),
                                                  StringIO(fasta3),
                                                  StringIO(fasta4)])
        self.assertEqual(result, '>hey\nagtcagtcagtc\n>you\nacctg')

    def testTwoFastaFilesNoCommons(self):
        fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta2 = '>all\nagtcagtcagtc\n>good\nacctg\n>over\natgggtc\n>here\n\
            atggctattgaactgtatct'

        result = fastamanipulations.fastaSubtract([StringIO(fasta1),
                                                  StringIO(fasta2)])
        self.assertEqual(result, '>hey\nagtcagtcagtc\n>you\nacctg')

    def testThreeFastaFilesNoCommons(self):
        fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta2 = '>all\nagtcagtcagtc\n>good\nacctg\n>over\natgggtc\n>here\n\
            atggctattgaactgtatct'
        fasta5 = '>this\nagtcagtcagtc\n>has\nacctg\n>many\natgggtc\n>rabbits\n\
            atggctattgaactgtatct'

        result = fastamanipulations.fastaSubtract([StringIO(fasta1),
                                                  StringIO(fasta2),
                                                  StringIO(fasta5)])
        self.assertEqual(result, None)
