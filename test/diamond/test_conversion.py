from six.moves import builtins
from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from ..mocking import mockOpen

from json import dumps

from dark.diamond.conversion import (JSONRecordsReader,
                                     DiamondTabularFormatReader)
from dark.reads import Reads, AARead


DIAMOND_RECORDS = """\
ACC94	INSV	29.6	0.003	1	EFII	178	295	SSSEV	175	285	295
ACC94	CASV	28.1	0.008	1	KLL	7	37	ITRV	9	39	300
ACC94	GoldenGate	28.1	0.009	1	IKSKL	7	35	EETSR	9	37	293
ACC94	GoldenGate	23.5	0.21	1	TIMSVV	177	240	DDMV	179	235	293
ACC94	InfluenzaC	25.0	0.084	1	LHVNYL	1	203	DEELKA	2	210	290
ACC94	InfluenzaC	18.5	9.1	1	SEIICEVLK	226	257	VETVAQ	20	45	290
ACC94	FERV	24.6	0.11	1	YSCFT-NSEK	176	276	LGKRMFC	152	243	270
AKAV	AKAV	634	0.0	1	GEPFSVYG	1	306	NIYGEP	1	306	306
AKAV	WYOV	401	7e-143	1	PFSVYGRF	1	306	GEPMS	1	294	294
BHAV	TAIV	28.1	0.008	1	PKELHGLI	14	118	SLKSKE	15	131	307
BHAV	SouthBay	28.1	0.009	1	CRPTF	4	293	EFVFIY	6	342	343
"""


class TestDiamondTabularFormatReader(TestCase):
    """
    Test the DiamondTabularFormatReader class.
    """

    def testDiamondParams(self):
        """
        When a DIAMOND file has been read, its parameters must be present
        in the reader instance.
        """
        mockOpener = mockOpen(read_data=DIAMOND_RECORDS)
        with patch.object(builtins, 'open', mockOpener):
            reader = DiamondTabularFormatReader('file.txt')
            list(reader.records())
            expectedRef = ('Buchfink et al., Fast and Sensitive Protein '
                           'Alignment using DIAMOND, Nature Methods, 12, '
                           '59–60 (2015)')
            self.assertEqual('DIAMOND', reader.params['application'])
            self.assertEqual('v0.8.23', reader.params['version'])
            self.assertEqual(expectedRef, reader.params['reference'])

    def testDiamondInput(self):
        """
        Test conversion of a chunk of DIAMOND output.
        """
        mockOpener = mockOpen(read_data=DIAMOND_RECORDS)
        with patch.object(builtins, 'open', mockOpener):
            reader = DiamondTabularFormatReader('file.txt')
            acc94, akav, bhav = list(reader.records())
            self.assertEqual(5, len(acc94['alignments']))
            self.assertEqual(2, len(akav['alignments']))
            self.assertEqual(2, len(bhav['alignments']))

_JSON_RECORDS = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59–60 '
                      '(2015)')
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 15225,
                'hsps': [
                    {
                        'sbjct_start': 5548,
                        'sbjct_end': 5450,
                        'bits': 165.393,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 2.73597e-40,
                        'query_end': 99,
                        'sbjct': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTCTAATGTGGTTG'),
                        'query': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTTTGTTGTGGTTG'),
                    }
                ]
            }
        ],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA'
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 15225,
                'hsps': [
                    {
                        'sbjct_start': 12320,
                        'sbjct_end': 12220,
                        'bits': 178.016,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACTATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 15191,
                'hsps': [
                    {
                        'sbjct_start': 12279,
                        'sbjct_end': 12179,
                        'bits': 123.915,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCTTGTGTAGATGAACCAACCCATGGTTTAGTGGGT'
                                  'CCTCTCTCACCACGTGTTAAACTGTTAACATTATATTTCTCT'
                                  'ATGATTATGCCACTAGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 15140,
                'hsps': [
                    {
                        'sbjct_start': 12213,
                        'sbjct_end': 12119,
                        'bits': 87.848,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('TTTTCT--TGGGTTGATGATCCTACCCATGGTTTAGTGGGAC'
                                  'CCCTCTCACCACGAGTCAAAAAATTGGAATTGTATTTTTCAA'
                                  'TTATTATACCACT'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA'
    },
    {
        'alignments': [],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA'
    },
    {
        'alignments': [],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA'
    }

]

_JSON_RECORDS_ONE_MIDDLE = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59–60 '
                      '(2015)')
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 15225,
                'hsps': [
                    {
                        'sbjct_start': 5548,
                        'sbjct_end': 5450,
                        'bits': 165.393,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 2.73597e-40,
                        'query_end': 99,
                        'sbjct': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTCTAATGTGGTTG'),
                        'query': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTTTGTTGTGGTTG'),
                    }
                ]
            }
        ],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA'
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 15225,
                'hsps': [
                    {
                        'sbjct_start': 12320,
                        'sbjct_end': 12220,
                        'bits': 178.016,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACTATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 15191,
                'hsps': [
                    {
                        'sbjct_start': 12279,
                        'sbjct_end': 12179,
                        'bits': 123.915,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCTTGTGTAGATGAACCAACCCATGGTTTAGTGGGT'
                                  'CCTCTCTCACCACGTGTTAAACTGTTAACATTATATTTCTCT'
                                  'ATGATTATGCCACTAGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 15140,
                'hsps': [
                    {
                        'sbjct_start': 12213,
                        'sbjct_end': 12119,
                        'bits': 87.848,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('TTTTCT--TGGGTTGATGATCCTACCCATGGTTTAGTGGGAC'
                                  'CCCTCTCACCACGAGTCAAAAAATTGGAATTGTATTTTTCAA'
                                  'TTATTATACCACT'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA'
    },
    {
        'alignments': [],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA'
    }

]

_JSON_RECORDS_ONE_END = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59–60 '
                      '(2015)')
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 15225,
                'hsps': [
                    {
                        'sbjct_start': 5548,
                        'sbjct_end': 5450,
                        'bits': 165.393,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 2.73597e-40,
                        'query_end': 99,
                        'sbjct': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTCTAATGTGGTTG'),
                        'query': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTTTGTTGTGGTTG'),
                    }
                ]
            }
        ],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA'
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 15225,
                'hsps': [
                    {
                        'sbjct_start': 12320,
                        'sbjct_end': 12220,
                        'bits': 178.016,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACTATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 15191,
                'hsps': [
                    {
                        'sbjct_start': 12279,
                        'sbjct_end': 12179,
                        'bits': 123.915,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCTTGTGTAGATGAACCAACCCATGGTTTAGTGGGT'
                                  'CCTCTCTCACCACGTGTTAAACTGTTAACATTATATTTCTCT'
                                  'ATGATTATGCCACTAGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 15140,
                'hsps': [
                    {
                        'sbjct_start': 12213,
                        'sbjct_end': 12119,
                        'bits': 87.848,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('TTTTCT--TGGGTTGATGATCCTACCCATGGTTTAGTGGGAC'
                                  'CCCTCTCACCACGAGTCAAAAAATTGGAATTGTATTTTTCAA'
                                  'TTATTATACCACT'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA'
    },
    {
        'alignments': [],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA'
    },

]

_JSON_RECORDS_ONE_START = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59–60 '
                      '(2015)')
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 15225,
                'hsps': [
                    {
                        'sbjct_start': 12320,
                        'sbjct_end': 12220,
                        'bits': 178.016,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACTATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 15191,
                'hsps': [
                    {
                        'sbjct_start': 12279,
                        'sbjct_end': 12179,
                        'bits': 123.915,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCTTGTGTAGATGAACCAACCCATGGTTTAGTGGGT'
                                  'CCTCTCTCACCACGTGTTAAACTGTTAACATTATATTTCTCT'
                                  'ATGATTATGCCACTAGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 15140,
                'hsps': [
                    {
                        'sbjct_start': 12213,
                        'sbjct_end': 12119,
                        'bits': 87.848,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('TTTTCT--TGGGTTGATGATCCTACCCATGGTTTAGTGGGAC'
                                  'CCCTCTCACCACGAGTCAAAAAATTGGAATTGTATTTTTCAA'
                                  'TTATTATACCACT'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA'
    },
    {
        'alignments': [],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA'
    },
    {
        'alignments': [],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA'
    }

]

_JSON_RECORDS_TWO_END = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59–60 '
                      '(2015)')
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 15225,
                'hsps': [
                    {
                        'sbjct_start': 5548,
                        'sbjct_end': 5450,
                        'bits': 165.393,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 2.73597e-40,
                        'query_end': 99,
                        'sbjct': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTCTAATGTGGTTG'),
                        'query': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTTTGTTGTGGTTG'),
                    }
                ]
            }
        ],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA'
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 15225,
                'hsps': [
                    {
                        'sbjct_start': 12320,
                        'sbjct_end': 12220,
                        'bits': 178.016,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACTATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 15191,
                'hsps': [
                    {
                        'sbjct_start': 12279,
                        'sbjct_end': 12179,
                        'bits': 123.915,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('TTTTTCTCTTGTGTAGATGAACCAACCCATGGTTTAGTGGGT'
                                  'CCTCTCTCACCACGTGTTAAACTGTTAACATTATATTTCTCT'
                                  'ATGATTATGCCACTAGC'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 15140,
                'hsps': [
                    {
                        'sbjct_start': 12213,
                        'sbjct_end': 12119,
                        'bits': 87.848,
                        'frame': [
                            1,
                            -1,
                        ],
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('TTTTCT--TGGGTTGATGATCCTACCCATGGTTTAGTGGGAC'
                                  'CCCTCTCACCACGAGTCAAAAAATTGGAATTGTATTTTTCAA'
                                  'TTATTATACCACT'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA'
    }

]

JSON = '\n'.join(dumps(record) for record in _JSON_RECORDS) + '\n'
JSON_ONE_MIDDLE = '\n'.join(dumps(record) for record in
                            _JSON_RECORDS_ONE_MIDDLE) + '\n'
JSON_ONE_END = '\n'.join(dumps(record) for record in
                         _JSON_RECORDS_ONE_END) + '\n'
JSON_ONE_START = '\n'.join(dumps(record) for record in
                           _JSON_RECORDS_ONE_START) + '\n'
JSON_TWO_END = '\n'.join(dumps(record) for record in
                         _JSON_RECORDS_TWO_END) + '\n'


class TestJSONRecordsReader(TestCase):
    """
    Test the JSONRecordsReader class.
    """
    def testCorrectNumberOfAlignments(self):
        """
        A JSONRecordsReader must return the expected number of alignments.
        """
        reads = Reads([
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA',
                'GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT'
                'AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC'),
        ])

        mockOpener = mockOpen(read_data=JSON)
        with patch.object(builtins, 'open', mockOpener):
            reader = JSONRecordsReader('file.json')
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(alignments))

    def testCorrectNumberOfAlignmentsMatchMissingMiddle(self):
        """
        A JSONRecordsReader must return the expected number of alignments, if
        a match is missing in the middle of the JSON file.
        """
        reads = Reads([
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA',
                'GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT'
                'AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC'),
        ])

        mockOpener = mockOpen(read_data=JSON_ONE_MIDDLE)
        with patch.object(builtins, 'open', mockOpener):
            reader = JSONRecordsReader('file.json')
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(alignments))

    def testCorrectNumberOfAlignmentsMatchMissingEnd(self):
        """
        A JSONRecordsReader must return the expected number of alignments, if
        the last read has no matches.
        """
        reads = Reads([
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA',
                'GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT'
                'AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC'),
        ])

        mockOpener = mockOpen(read_data=JSON_ONE_END)
        with patch.object(builtins, 'open', mockOpener):
            reader = JSONRecordsReader('file.json')
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(alignments))

    def testCorrectNumberOfAlignmentsMatchMissingStart(self):
        """
        A JSONRecordsReader must return the expected number of alignments, if
        the first read has no matches.
        """
        reads = Reads([
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA',
                'GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT'
                'AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC'),
        ])

        mockOpener = mockOpen(read_data=JSON_ONE_START)
        with patch.object(builtins, 'open', mockOpener):
            reader = JSONRecordsReader('file.json')
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(alignments))

    def testCorrectNumberOfAlignmentsTwoMatchesMissingEnd(self):
        """
        A JSONRecordsReader must return the expected number of alignments, if
        two reads at the end don't have any matches.
        """
        reads = Reads([
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234 1:N:0:TGACCA',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287 1:N:0:TGACCA',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488 1:N:0:TGACCA',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512 1:N:0:TGACCA',
                'GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT'
                'AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC'),
        ])

        mockOpener = mockOpen(read_data=JSON_TWO_END)
        with patch.object(builtins, 'open', mockOpener):
            reader = JSONRecordsReader('file.json')
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(alignments))

    def testCorrectNumberOfAlignmentsWhenReadIdsAreAbbreviated(self):
        """
        A JSONRecordsReader must return the expected number of alignments
        when read ids are truncated at the first space. That is, the DIAMOND
        output has query names that are long and which contain a space but
        the reads in the FASTA have just the first part of those names (up to
        the first space).
        """
        reads = Reads([
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:9489:4234',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:19964:6287',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:11488:7488',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'BIOMICS-HISEQTP:140:HJFH5BCXX:1:1101:14734:7512',
                'GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT'
                'AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC'),
        ])

        mockOpener = mockOpen(read_data=JSON)
        with patch.object(builtins, 'open', mockOpener):
            reader = JSONRecordsReader('file.json')
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(alignments))
