from unittest import TestCase
from cStringIO import StringIO
from Bio import SeqIO
from dark import subtract


class TestFastaSubtract(TestCase):
    """
    Test the fastaSubtract() function.
    """

    def testTwoFastaFilesInvertFalse(self):
        fasta1 = '>hey\nagtcagtcagtc\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta3 = '>hey\nagtcagtcagtc\n>this\natgggtc\n>is\n\
            atggctattgaactgtatct'
        testResult = []
        testResult.append(SeqIO.read(StringIO('>hey\nagtcagtcagtc'), 'fasta'))

        result = subtract.fastaSubtract([StringIO(fasta1),
                                        StringIO(fasta3)], invert=False)

        self.assertEqual(result, testResult)

    def testTwoFastaFilesInvertTrue(self):
        fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta3 = '>hey\nagtcagtcagtc\n>you\nacctg\n>this\natgggtc\n>is\n\
            atggctattgaactgtatct'
        testResult = []
        testResult.append(SeqIO.read(StringIO('>how\natgggtc'), 'fasta'))
        testResult.append(SeqIO.read(StringIO('>are\natggctattgaactgtatct'),
                                     'fasta'))
        testResult.append(SeqIO.read(StringIO('>this\natgggtc'), 'fasta'))
        testResult.append(SeqIO.read(StringIO('>is\natggctattgaactgtatct'),
                                     'fasta'))

        result = subtract.fastaSubtract([StringIO(fasta1),
                                        StringIO(fasta3)], invert=True)
        self.assertEqual(result, testResult)

    def testThreeFastaFilesInvertFalse(self):
        fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta3 = '>hey\nagtcagtcagtc\n>you\nacctg\n>this\natgggtc\n>is\n\
            atggctattgaactgtatct'
        fasta4 = '>hey\nagtcagtcagtc\n>you\nacctg\n>nothing\natgggtc\n>there\n\
            atggctattgaactgtatct'

        testResult = []
        testResult.append(SeqIO.read(StringIO('>hey\nagtcagtcagtc'), 'fasta'))
        testResult.append(SeqIO.read(StringIO('>you\nacctg'), 'fasta'))

        result = subtract.fastaSubtract([StringIO(fasta1),
                                        StringIO(fasta3),
                                        StringIO(fasta4)], invert=False)
        self.assertEqual(result, testResult)

    def testThreeFastaFilesInvertTrue(self):
        fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta3 = '>hey\nagtcagtcagtc\n>you\nacctg\n>this\natgggtc\n>is\n\
            atggctattgaactgtatct'
        fasta4 = '>hey\nagtcagtcagtc\n>you\nacctg\n>nothing\natgggtc\n>there\n\
            atggctattgaactgtatct'
        testResult = []
        testResult.append(SeqIO.read(StringIO('>how\natgggtc'), 'fasta'))
        testResult.append(SeqIO.read(StringIO('>are\natggctattgaactgtatct'),
                                     'fasta'))
        testResult.append(SeqIO.read(StringIO('>this\natgggtc'), 'fasta'))
        testResult.append(SeqIO.read(StringIO('>is\natggctattgaactgtatct'),
                                     'fasta'))
        testResult.append(SeqIO.read(StringIO('>nothing\natgggtc'), 'fasta'))
        testResult.append(SeqIO.read(StringIO('>there\natggctattgaactgtatct'),
                                     'fasta'))

        result = subtract.fastaSubtract([StringIO(fasta1),
                                        StringIO(fasta3),
                                        StringIO(fasta4)], invert=True)
        self.assertEqual(result, testResult)

    def testTwoFastaFilesNoCommons(self):
        fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta2 = '>all\nagtcagtcagtc\n>good\nacctg\n>over\natgggtc\n>here\n\
            atggctattgaactgtatct'

        result = subtract.fastaSubtract([StringIO(fasta1),
                                        StringIO(fasta2)])
        self.assertEqual(result, [])

    def testThreeFastaFilesNoCommons(self):
        fasta1 = '>hey\nagtcagtcagtc\n>you\nacctg\n>how\natgggtc\n>are\n\
            atggctattgaactgtatct'
        fasta2 = '>all\nagtcagtcagtc\n>good\nacctg\n>over\natgggtc\n>here\n\
            atggctattgaactgtatct'
        fasta5 = '>this\nagtcagtcagtc\n>has\nacctg\n>many\natgggtc\n>rabbits\n\
            atggctattgaactgtatct'

        result = subtract.fastaSubtract([StringIO(fasta1),
                                        StringIO(fasta2),
                                        StringIO(fasta5)])
        self.assertEqual(result, [])
