from unittest import TestCase
from Bio import SeqIO
from dark.smith_waterman import smith_waterman
from cStringIO import StringIO

class TestSmithWaterman(TestCase):
    """
    Test the smith_waterman function.
    """

    def testZeroInput(self):
    	seq1 = ''
    	seq2 = ''
    	result = smith_waterman(seq1,seq2)
    	self.assertEqual(result, [])

    def testOneInput(self):
    	seq1 = 'agtcagtcagtc'
    	seq2 = ''
    	result = smith_waterman(seq1,seq2)
    	self.assertEqual(result, [])

    def testTwoInput(self):
    	seq1 = 'agtcagtcagtc'
    	seq2 = 'cgaatcg'
    	result = smith_waterman(seq1,seq2)
    	self.assertEqual(result, ['agtcagtcagtc','      || |  ','--cgaatc-g--'])

