from unittest import TestCase

from dark.blast.score import bitScoreToEValue, eValueToBitScore


class TestBitScoreToEValue(TestCase):
    """
    Tests for the bitScoreToEValue function.
    """

    def testTrivial(self):
        """
        For a bit score of 1.0 on a database of size 1 with a query length
        1, and no length adjustment, we should get an e-value of 0.5
        """
        self.assertEqual(
            0.5,
            bitScoreToEValue(
                bitScore=1.0,
                dbSize=1,
                dbSequenceCount=1,
                queryLength=1,
                lengthAdjustment=0,
            ),
        )

    def testErwiniaPhage(self):
        """
        Test values corresponding to the following observed BLAST match against
        Erwinia phage phiEaH2, complete genome.

        <Iteration>
          <Iteration_iter-num>3</Iteration_iter-num>
          <Iteration_query-ID>Query_3</Iteration_query-ID>
          <Iteration_query-def>SK7F6:834:2401</Iteration_query-def>
          <Iteration_query-len>111</Iteration_query-len>
        <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|431810571|ref|NC_019929.1|</Hit_id>
          <Hit_def>Erwinia phage phiEaH2, complete genome</Hit_def>
          <Hit_accession>NC_019929</Hit_accession>
          <Hit_len>243050</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>37.3537</Hsp_bit-score>
              <Hsp_score>40</Hsp_score>
              <Hsp_evalue>0.0813089</Hsp_evalue>
              <Hsp_query-from>86</Hsp_query-from>
              <Hsp_query-to>108</Hsp_query-to>
              <Hsp_hit-from>33712</Hsp_hit-from>
              <Hsp_hit-to>33690</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>-1</Hsp_hit-frame>
              <Hsp_identity>22</Hsp_identity>
              <Hsp_positive>22</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>23</Hsp_align-len>
              <Hsp_qseq>GATAACTTCCAGGAATTCGTCCA</Hsp_qseq>
              <Hsp_hseq>GATAACTTCCAGGAATTCGTTCA</Hsp_hseq>
              <Hsp_midline>|||||||||||||||||||| ||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        </Iteration_hits>
          <Iteration_stat>
            <Statistics>
              <Statistics_db-num>5660</Statistics_db-num>
              <Statistics_db-len>168142520</Statistics_db-len>
              <Statistics_hsp-len>26</Statistics_hsp-len>
              <Statistics_eff-space>14279605600</Statistics_eff-space>
              <Statistics_kappa>0.41</Statistics_kappa>
              <Statistics_lambda>0.625</Statistics_lambda>
              <Statistics_entropy>0.78</Statistics_entropy>
            </Statistics>
          </Iteration_stat>
        </Iteration>
        """
        self.assertAlmostEqual(
            0.0813089,
            bitScoreToEValue(
                bitScore=37.3537,
                dbSize=168142520,
                dbSequenceCount=5660,
                queryLength=111,
                lengthAdjustment=26,
            ),
            places=4,
        )

    def testParameciumBursaria(self):
        """
        Test values corresponding to the following observed BLAST match against
        Paramecium bursaria Chlorella virus NE-JV-1, partial genome.
        (Long lines are continued with a backslash to keep Python linters
        quiet.)

        <Iteration>
          <Iteration_iter-num>14</Iteration_iter-num>
          <Iteration_query-ID>Query_14</Iteration_query-ID>
          <Iteration_query-def>SK7F6:99:509</Iteration_query-def>
          <Iteration_query-len>172</Iteration_query-len>
        <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|448932590|gb|JX997176.1|</Hit_id>
          <Hit_def>Paramecium bursaria Chlorella virus NE-JV-1, partial \
              genome</Hit_def>
          <Hit_accession>JX997176</Hit_accession>
          <Hit_len>326559</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>42.7638</Hsp_bit-score>
              <Hsp_score>46</Hsp_score>
              <Hsp_evalue>0.0359052</Hsp_evalue>
              <Hsp_query-from>93</Hsp_query-from>
              <Hsp_query-to>164</Hsp_query-to>
              <Hsp_hit-from>140106</Hsp_hit-from>
              <Hsp_hit-to>140038</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>-1</Hsp_hit-frame>
              <Hsp_identity>55</Hsp_identity>
              <Hsp_positive>55</Hsp_positive>
              <Hsp_gaps>3</Hsp_gaps>
              <Hsp_align-len>72</Hsp_align-len>
              <Hsp_qseq>AGGAATCCCTAAATGAAGTCCAAGAAGAAATCCCCTGAAGGAATTATAAAGAGA\
                  AATCCCTGAAGTAATTCC</Hsp_qseq>
              <Hsp_hseq>AGAAATCCCTGAA-GAAATCCCTGAAGACATCCC-TGAAGAAATCCCTGA-AGA\
                  AATCCCTGAAGAAATCCC</Hsp_hseq>
              <Hsp_midline>|| ||||||| || ||| |||  ||||| ||||| ||||| |||     | \
                  |||||||||||||| ||| ||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        </Iteration_hits>
          <Iteration_stat>
            <Statistics>
              <Statistics_db-num>1456080</Statistics_db-num>
              <Statistics_db-len>1931895878</Statistics_db-len>
              <Statistics_hsp-len>30</Statistics_hsp-len>
              <Statistics_eff-space>268126313876</Statistics_eff-space>
              <Statistics_kappa>0.41</Statistics_kappa>
              <Statistics_lambda>0.625</Statistics_lambda>
              <Statistics_entropy>0.78</Statistics_entropy>
            </Statistics>
          </Iteration_stat>
        </Iteration>
        """
        self.assertAlmostEqual(
            0.0359052,
            bitScoreToEValue(
                bitScore=42.7638,
                dbSize=1931895878,
                dbSequenceCount=1456080,
                queryLength=172,
                lengthAdjustment=30,
            ),
            places=4,
        )


class TestEValueToBitScore(TestCase):
    """
    Tests for the eValueToBitScore function.
    """

    def testTrivial(self):
        """
        For an e-value of 1.0 on a database of size 1 with total query length
        1, we should get a bit score of 0.0
        """
        self.assertEqual(
            0.0,
            eValueToBitScore(
                eValue=1.0,
                dbSize=1,
                dbSequenceCount=1,
                queryLength=1,
                lengthAdjustment=0,
            ),
        )

    def testErwiniaPhage(self):
        """
        Test values corresponding to the following observed BLAST match against
        Erwinia phage phiEaH2, complete genome.

        <Iteration>
          <Iteration_iter-num>3</Iteration_iter-num>
          <Iteration_query-ID>Query_3</Iteration_query-ID>
          <Iteration_query-def>SK7F6:834:2401</Iteration_query-def>
          <Iteration_query-len>111</Iteration_query-len>
        <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|431810571|ref|NC_019929.1|</Hit_id>
          <Hit_def>Erwinia phage phiEaH2, complete genome</Hit_def>
          <Hit_accession>NC_019929</Hit_accession>
          <Hit_len>243050</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>37.3537</Hsp_bit-score>
              <Hsp_score>40</Hsp_score>
              <Hsp_evalue>0.0813089</Hsp_evalue>
              <Hsp_query-from>86</Hsp_query-from>
              <Hsp_query-to>108</Hsp_query-to>
              <Hsp_hit-from>33712</Hsp_hit-from>
              <Hsp_hit-to>33690</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>-1</Hsp_hit-frame>
              <Hsp_identity>22</Hsp_identity>
              <Hsp_positive>22</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>23</Hsp_align-len>
              <Hsp_qseq>GATAACTTCCAGGAATTCGTCCA</Hsp_qseq>
              <Hsp_hseq>GATAACTTCCAGGAATTCGTTCA</Hsp_hseq>
              <Hsp_midline>|||||||||||||||||||| ||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        </Iteration_hits>
          <Iteration_stat>
            <Statistics>
              <Statistics_db-num>5660</Statistics_db-num>
              <Statistics_db-len>168142520</Statistics_db-len>
              <Statistics_hsp-len>26</Statistics_hsp-len>
              <Statistics_eff-space>14279605600</Statistics_eff-space>
              <Statistics_kappa>0.41</Statistics_kappa>
              <Statistics_lambda>0.625</Statistics_lambda>
              <Statistics_entropy>0.78</Statistics_entropy>
            </Statistics>
          </Iteration_stat>
        </Iteration>
        """
        self.assertAlmostEqual(
            37.3537,
            eValueToBitScore(
                eValue=0.0813089,
                dbSize=168142520,
                dbSequenceCount=5660,
                queryLength=111,
                lengthAdjustment=26,
            ),
            places=4,
        )

    def testParameciumBursaria(self):
        """
        Test values corresponding to the following observed BLAST match against
        Paramecium bursaria Chlorella virus NE-JV-1, partial genome.
        (Long lines are continued with a backslash to keep Python linters
        quiet.)

        <Iteration>
          <Iteration_iter-num>14</Iteration_iter-num>
          <Iteration_query-ID>Query_14</Iteration_query-ID>
          <Iteration_query-def>SK7F6:99:509</Iteration_query-def>
          <Iteration_query-len>172</Iteration_query-len>
        <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|448932590|gb|JX997176.1|</Hit_id>
          <Hit_def>Paramecium bursaria Chlorella virus NE-JV-1, partial \
              genome</Hit_def>
          <Hit_accession>JX997176</Hit_accession>
          <Hit_len>326559</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>42.7638</Hsp_bit-score>
              <Hsp_score>46</Hsp_score>
              <Hsp_evalue>0.0359052</Hsp_evalue>
              <Hsp_query-from>93</Hsp_query-from>
              <Hsp_query-to>164</Hsp_query-to>
              <Hsp_hit-from>140106</Hsp_hit-from>
              <Hsp_hit-to>140038</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>-1</Hsp_hit-frame>
              <Hsp_identity>55</Hsp_identity>
              <Hsp_positive>55</Hsp_positive>
              <Hsp_gaps>3</Hsp_gaps>
              <Hsp_align-len>72</Hsp_align-len>
              <Hsp_qseq>AGGAATCCCTAAATGAAGTCCAAGAAGAAATCCCCTGAAGGAATTATAAAGAGA\
                  AATCCCTGAAGTAATTCC</Hsp_qseq>
              <Hsp_hseq>AGAAATCCCTGAA-GAAATCCCTGAAGACATCCC-TGAAGAAATCCCTGA-AGA\
                  AATCCCTGAAGAAATCCC</Hsp_hseq>
              <Hsp_midline>|| ||||||| || ||| |||  ||||| ||||| ||||| |||     | \
                  |||||||||||||| ||| ||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        </Iteration_hits>
          <Iteration_stat>
            <Statistics>
              <Statistics_db-num>1456080</Statistics_db-num>
              <Statistics_db-len>1931895878</Statistics_db-len>
              <Statistics_hsp-len>30</Statistics_hsp-len>
              <Statistics_eff-space>268126313876</Statistics_eff-space>
              <Statistics_kappa>0.41</Statistics_kappa>
              <Statistics_lambda>0.625</Statistics_lambda>
              <Statistics_entropy>0.78</Statistics_entropy>
            </Statistics>
          </Iteration_stat>
        </Iteration>
        """
        self.assertAlmostEqual(
            42.7638,
            eValueToBitScore(
                eValue=0.0359052,
                dbSize=1931895878,
                dbSequenceCount=1456080,
                queryLength=172,
                lengthAdjustment=30,
            ),
            places=4,
        )
