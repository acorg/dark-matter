import builtins
from json import dumps
from unittest import TestCase
from unittest.mock import mock_open, patch

from dark.blast.conversion import JSONRecordsReader, XMLRecordsReader
from dark.reads import DNARead, Reads

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
              <Hsp_positive>17</Hsp_positive>
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
              <Hsp_identity>18</Hsp_identity>
              <Hsp_positive>19</Hsp_positive>
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
        data = RECORD
        with patch.object(builtins, "open", mock_open(read_data=data)):
            reader = XMLRecordsReader("file.xml")
            list(reader.records())
            self.assertEqual("BLASTN", reader.params["application"])
            self.assertEqual("virus-nt-20130719", reader.params["database"])

    def testXMLInput(self):
        """
        Test conversion of a chunk of BLAST XML. This is highly incomplete
        in what it tests.
        """
        mockOpener = mock_open(read_data=RECORD)
        with patch.object(builtins, "open", mockOpener):
            reader = XMLRecordsReader("file.xml")
            record1, record2 = list(reader.records())
            self.assertEqual(0, len(record1.alignments))
            self.assertEqual(2, len(record2.alignments))


_JSON_RECORDS = [
    {
        "application": "BLASTN",
        "effective_search_space": 16800583950.0,
        "num_sequences_in_database": 7194,
        "version": "2.2.30+",
        "query_letters": 101,
        "database_letters": None,
        "gap_penalties": [
            5,
            2,
        ],
        "database": "virus-refseq-20160316",
        "database_sequences": 7194,
        "query": "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA",
        "hsps_prelim_gapped_attemped": None,
        "sc_mismatch": -3,
        "effective_hsp_length": 26,
        "threshold": None,
        "num_seqs_better_e": None,
        "num_hits": None,
        "hsps_no_gap": None,
        "hsps_gapped": None,
        "matrix": "",
        "database_length": 224194830,
        "sc_match": 2,
        "gap_trigger": [
            None,
            None,
        ],
        "num_good_extends": None,
        "ka_params": [
            0.625,
            0.41,
            0.78,
        ],
        "gap_x_dropoff_final": [
            None,
            None,
        ],
        "effective_search_space_used": None,
        "date": "",
        "query_id": "Query_1",
        "database_name": [],
        "num_sequences": None,
        "effective_database_length": None,
        "reference": "Stephen F. Altschul, ...",
        "dropoff_1st_pass": [
            None,
            None,
        ],
        "hsps_prelim_gapped": None,
        "frameshift": [
            None,
            None,
        ],
        "window_size": None,
        "query_length": 101,
        "num_letters_in_database": 224194830,
        "gapped": 0,
        "gap_x_dropoff": [
            None,
            None,
        ],
        "blast_cutoff": [
            None,
            None,
        ],
        "posted_date": [],
        "effective_query_length": None,
        "ka_params_gap": [
            None,
            None,
            None,
        ],
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 15225,
                "hsps": [
                    {
                        "sbjct_start": 5548,
                        "sbjct_end": 5450,
                        "bits": 165.393,
                        "frame": [
                            1,
                            -1,
                        ],
                        "identicalCount": 38,
                        "positiveCount": 77,
                        "query_start": 1,
                        "expect": 2.73597e-40,
                        "query_end": 99,
                        "sbjct": (
                            "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT"
                            "TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT"
                            "GTTCTAATGTGGTTG"
                        ),
                        "query": (
                            "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT"
                            "TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT"
                            "GTTTTGTTGTGGTTG"
                        ),
                    }
                ],
            }
        ],
        "query": "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 15225,
                "hsps": [
                    {
                        "sbjct_start": 12320,
                        "sbjct_end": 12220,
                        "bits": 178.016,
                        "frame": [
                            1,
                            -1,
                        ],
                        "identicalCount": 380,
                        "positiveCount": 770,
                        "query_start": 1,
                        "expect": 4.33545e-44,
                        "query_end": 101,
                        "sbjct": (
                            "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT"
                            "CCTCTTTCACCACGAGTTAAACTATTAACATTATATTTTTCT"
                            "ATAATTATACCACTGGC"
                        ),
                        "query": (
                            "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT"
                            "CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT"
                            "ATAATTATACCACTGGC"
                        ),
                    }
                ],
            },
            {
                "title": "gi|9629367|ref|NC_001803.1| RSV",
                "length": 15191,
                "hsps": [
                    {
                        "sbjct_start": 12279,
                        "sbjct_end": 12179,
                        "bits": 123.915,
                        "frame": [
                            1,
                            -1,
                        ],
                        "identicalCount": 3,
                        "positiveCount": 7,
                        "query_start": 1,
                        "expect": 8.37678e-28,
                        "query_end": 101,
                        "sbjct": (
                            "TTTTTCTCTTGTGTAGATGAACCAACCCATGGTTTAGTGGGT"
                            "CCTCTCTCACCACGTGTTAAACTGTTAACATTATATTTCTCT"
                            "ATGATTATGCCACTAGC"
                        ),
                        "query": (
                            "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT"
                            "CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT"
                            "ATAATTATACCACTGGC"
                        ),
                    }
                ],
            },
            {
                "title": "gi|9631267|ref|NC_001989.1| Bovine RSV",
                "length": 15140,
                "hsps": [
                    {
                        "sbjct_start": 12213,
                        "sbjct_end": 12119,
                        "bits": 87.848,
                        "frame": [
                            1,
                            -1,
                        ],
                        "identicalCount": 3800,
                        "positiveCount": 7700,
                        "query_start": 2,
                        "expect": 6.03169e-17,
                        "query_end": 98,
                        "sbjct": (
                            "TTTTCT--TGGGTTGATGATCCTACCCATGGTTTAGTGGGAC"
                            "CCCTCTCACCACGAGTCAAAAAATTGGAATTGTATTTTTCAA"
                            "TTATTATACCACT"
                        ),
                        "query": (
                            "TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC"
                            "CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA"
                            "TAATTATACCACT"
                        ),
                    }
                ],
            },
        ],
        "query": "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA",
    },
    {
        "alignments": [],
        "query": "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA",
    },
    {
        "alignments": [],
        "query": "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA",
    },
]

JSON = "\n".join(dumps(record) for record in _JSON_RECORDS) + "\n"


class TestJSONRecordsReader(TestCase):
    """
    Test the JSONRecordsReader class.
    """

    READS = Reads(
        [
            DNARead(
                "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA",
                "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG"
                "GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG",
            ),
            DNARead(
                "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA",
                "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC"
                "ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC",
            ),
            DNARead(
                "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA",
                "ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA"
                "GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT",
            ),
            DNARead(
                "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA",
                "GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT"
                "AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC",
            ),
        ]
    )

    # TODO: This class is quite incomplete. It was added long after the
    # original code was written, due to laziness at the time. Additional
    # aspects of the returned alignments should be tested.

    def testCorrectNumberOfReadAlignments(self):
        """
        A JSONRecordsReader must return the expected number of read alignments.
        """
        mockOpener = mock_open(read_data=JSON)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            alignments = list(reader.readAlignments(self.READS))
            self.assertEqual(4, len(alignments))

    def testCorrectNumberOfAlignmentsWhenReadIdsAreAbbreviated(self):
        """
        A JSONRecordsReader must return the expected number of alignments
        when read ids are truncated at the first space. That is, the BLAST
        output has query names that are long and which contain a space but
        the reads in the FASTA have just the first part of those names (up to
        the first space).
        """
        reads = Reads(
            [
                DNARead(
                    "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234",
                    "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG"
                    "GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG",
                ),
                DNARead(
                    "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287",
                    "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC"
                    "ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC",
                ),
                DNARead(
                    "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488",
                    "ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA"
                    "GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT",
                ),
                DNARead(
                    "BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512",
                    "GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT"
                    "AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC",
                ),
            ]
        )

        mockOpener = mock_open(read_data=JSON)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            readAlignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(readAlignments))

    def testIdentity(self):
        """
        The identity value must be read correctly from the HSPs.
        """
        mockOpener = mock_open(read_data=JSON)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            readAlignments = list(reader.readAlignments(self.READS))
            self.assertEqual(38, readAlignments[0][0].hsps[0].identicalCount)
            self.assertEqual(77, readAlignments[0][0].hsps[0].positiveCount)
            self.assertEqual(380, readAlignments[1][0].hsps[0].identicalCount)
            self.assertEqual(770, readAlignments[1][0].hsps[0].positiveCount)
            self.assertEqual(3, readAlignments[1][1].hsps[0].identicalCount)
            self.assertEqual(7, readAlignments[1][1].hsps[0].positiveCount)
            self.assertEqual(3800, readAlignments[1][2].hsps[0].identicalCount)
            self.assertEqual(7700, readAlignments[1][2].hsps[0].positiveCount)
