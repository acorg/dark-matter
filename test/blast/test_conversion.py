from six.moves import builtins
from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from ..mocking import mockOpen

from dark.blast.conversion import XMLRecordsReader

RECORD = """\
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" \
"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.2.27+</BlastOutput_version>
  <BlastOutput_reference>[Snip]</BlastOutput_reference>
  <BlastOutput_db>virus-nt-20130719</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>ICUR3MX01AGWKS</BlastOutput_query-def>
  <BlastOutput_query-len>63</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>10</Parameters_expect>
      <Parameters_sc-match>2</Parameters_sc-match>
      <Parameters_sc-mismatch>-3</Parameters_sc-mismatch>
      <Parameters_gap-open>5</Parameters_gap-open>
      <Parameters_gap-extend>2</Parameters_gap-extend>
      <Parameters_filter>L;m;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>ICUR3MX01AGWKS</Iteration_query-def>
      <Iteration_query-len>63</Iteration_query-len>
      <Iteration_hits></Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1331289</Statistics_db-num>
          <Statistics_db-len>1634236772</Statistics_db-len>
          <Statistics_hsp-len>28</Statistics_hsp-len>
          <Statistics_eff-space>55893623800</Statistics_eff-space>
          <Statistics_kappa>0.41</Statistics_kappa>
          <Statistics_lambda>0.625</Statistics_lambda>
          <Statistics_entropy>0.78</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
      <Iteration_message>No hits found</Iteration_message>
    </Iteration>
    <Iteration>
      <Iteration_iter-num>2</Iteration_iter-num>
      <Iteration_query-ID>Query_2</Iteration_query-ID>
      <Iteration_query-def>ICUR3MX01C58U5</Iteration_query-def>
      <Iteration_query-len>21</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>gi|157694982|gb|EF990688.1|</Hit_id>
          <Hit_def>Rotavirus A strain RVA/Pig-tc/VEN/A131/1988</Hit_def>
          <Hit_accession>EF990688</Hit_accession>
          <Hit_len>1059</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>30.1402</Hsp_bit-score>
              <Hsp_score>32</Hsp_score>
              <Hsp_evalue>4.0824</Hsp_evalue>
              <Hsp_query-from>6</Hsp_query-from>
              <Hsp_query-to>21</Hsp_query-to>
              <Hsp_hit-from>738</Hsp_hit-from>
              <Hsp_hit-to>753</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>16</Hsp_identity>
              <Hsp_positive>16</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>16</Hsp_align-len>
              <Hsp_qseq>ATTCATCAGTAGCAAT</Hsp_qseq>
              <Hsp_hseq>ATTCATCAGTAGCAAT</Hsp_hseq>
              <Hsp_midline>||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>gi|157694972|gb|EF990692.1|</Hit_id>
          <Hit_def>Rotavirus A strain RVA/Pig-tc/VEN/A411/1989</Hit_def>
          <Hit_accession>EF990692</Hit_accession>
          <Hit_len>1059</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>30.1402</Hsp_bit-score>
              <Hsp_score>32</Hsp_score>
              <Hsp_evalue>4.0824</Hsp_evalue>
              <Hsp_query-from>6</Hsp_query-from>
              <Hsp_query-to>21</Hsp_query-to>
              <Hsp_hit-from>738</Hsp_hit-from>
              <Hsp_hit-to>753</Hsp_hit-to>
              <Hsp_query-frame>1</Hsp_query-frame>
              <Hsp_hit-frame>1</Hsp_hit-frame>
              <Hsp_identity>16</Hsp_identity>
              <Hsp_positive>16</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>16</Hsp_align-len>
              <Hsp_qseq>ATTCATCAGTAGCAAT</Hsp_qseq>
              <Hsp_hseq>ATTCATCAGTAGCAAT</Hsp_hseq>
              <Hsp_midline>||||||||||||||||</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1331289</Statistics_db-num>
          <Statistics_db-len>1634236772</Statistics_db-len>
          <Statistics_hsp-len>18</Statistics_hsp-len>
          <Statistics_eff-space>4830820710</Statistics_eff-space>
          <Statistics_kappa>0.41</Statistics_kappa>
          <Statistics_lambda>0.625</Statistics_lambda>
          <Statistics_entropy>0.78</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


class TestXMLRecordsReader(TestCase):
    """
    Test the XMLRecordsReader class.
    """

    def testParams(self):
        """
        When a BLAST XML file has been read, its parameters must be present
        in the reader instance. We only test a subset of the parameters.
        """
        mockOpener = mockOpen(read_data=RECORD)
        with patch.object(builtins, 'open', mockOpener):
            reader = XMLRecordsReader('file.xml')
            list(reader.records())
            self.assertEqual('BLASTN', reader.params['application'])
            self.assertEqual('virus-nt-20130719', reader.params['database'])

    def testXMLInput(self):
        """
        Test conversion of a chunk of BLAST XML. This is highly incomplete
        in what it tests.
        """
        mockOpener = mockOpen(read_data=RECORD)
        with patch.object(builtins, 'open', mockOpener):
            reader = XMLRecordsReader('file.xml')
            record1, record2 = list(reader.records())
            self.assertEqual(0, len(record1.alignments))
            self.assertEqual(2, len(record2.alignments))
