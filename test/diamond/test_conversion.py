from six.moves import builtins
from unittest import TestCase
from io import BytesIO, StringIO
import bz2file
from bz2 import compress

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

DIAMOND_RECORDS_DUMPED = '\n'.join([
    dumps({
        "reference": ("Buchfink et al., Fast and Sensitive "
                      "Protein Alignment using DIAMOND, Nature Methods, "
                      "12, 59-60 (2015)"),
        "task": "blastx",
        "version": "v0.8.23"
    }, sort_keys=True),
    dumps({
        "alignments": [
            {
                "hsps": [
                    {
                        "bits": 29.6,
                        "expect": 0.003,
                        "frame": 1,
                        "query": "EFII",
                        "query_end": 295,
                        "query_start": 178,
                        "sbjct": "SSSEV",
                        "sbjct_end": 285,
                        "sbjct_start": 175
                    }
                ],
                "length": 295,
                "title": "INSV"
            },
            {
                "hsps": [
                    {
                        "bits": 28.1,
                        "expect": 0.008,
                        "frame": 1,
                        "query": "KLL",
                        "query_end": 37,
                        "query_start": 7,
                        "sbjct": "ITRV",
                        "sbjct_end": 39,
                        "sbjct_start": 9
                    }
                ],
                "length": 300,
                "title": "CASV"
            },
            {
                "hsps": [
                    {
                        "bits": 28.1,
                        "expect": 0.009,
                        "frame": 1,
                        "query": "IKSKL",
                        "query_end": 35,
                        "query_start": 7,
                        "sbjct": "EETSR",
                        "sbjct_end": 37,
                        "sbjct_start": 9
                    },
                    {
                        "bits": 23.5,
                        "expect": 0.21,
                        "frame": 1,
                        "query": "TIMSVV",
                        "query_end": 240,
                        "query_start": 177,
                        "sbjct": "DDMV",
                        "sbjct_end": 235,
                        "sbjct_start": 179
                    }
                ],
                "length": 293,
                "title": "GoldenGate"
            },
            {
                "hsps": [
                    {
                        "bits": 25.0,
                        "expect": 0.084,
                        "frame": 1,
                        "query": "LHVNYL",
                        "query_end": 203,
                        "query_start": 1,
                        "sbjct": "DEELKA",
                        "sbjct_end": 210,
                        "sbjct_start": 2
                    },
                    {
                        "bits": 18.5,
                        "expect": 9.1,
                        "frame": 1,
                        "query": "SEIICEVLK",
                        "query_end": 257,
                        "query_start": 226,
                        "sbjct": "VETVAQ",
                        "sbjct_end": 45,
                        "sbjct_start": 20
                    }
                ],
                "length": 290,
                "title": "InfluenzaC"
            },
            {
                "hsps": [
                    {
                        "bits": 24.6,
                        "expect": 0.11,
                        "frame": 1,
                        "query": "YSCFT-NSEK",
                        "query_end": 276,
                        "query_start": 176,
                        "sbjct": "LGKRMFC",
                        "sbjct_end": 243,
                        "sbjct_start": 152
                    }
                ],
                "length": 270,
                "title": "FERV"
            }
        ],
        "query": "ACC94"
    }, sort_keys=True),
    dumps({
        "alignments": [
            {
                "hsps": [
                    {
                        "bits": 634.0,
                        "expect": 0.0,
                        "frame": 1,
                        "query": "GEPFSVYG",
                        "query_end": 306,
                        "query_start": 1,
                        "sbjct": "NIYGEP",
                        "sbjct_end": 306,
                        "sbjct_start": 1
                    }
                ],
                "length": 306,
                "title": "AKAV"
            },
            {
                "hsps": [
                    {
                        "bits": 401.0,
                        "expect": 7e-143,
                        "frame": 1,
                        "query": "PFSVYGRF",
                        "query_end": 306,
                        "query_start": 1,
                        "sbjct": "GEPMS",
                        "sbjct_end": 294,
                        "sbjct_start": 1
                    }
                ],
                "length": 294,
                "title": "WYOV"
            }
        ],
        "query": "AKAV"
    }, sort_keys=True),
    dumps({
        "alignments": [
            {
                "hsps": [
                    {
                        "bits": 28.1,
                        "expect": 0.008,
                        "frame": 1,
                        "query": "PKELHGLI",
                        "query_end": 118,
                        "query_start": 14,
                        "sbjct": "SLKSKE",
                        "sbjct_end": 131,
                        "sbjct_start": 15
                    }
                ],
                "length": 307,
                "title": "TAIV"
            },
            {
                "hsps": [
                    {
                        "bits": 28.1,
                        "expect": 0.009,
                        "frame": 1,
                        "query": "CRPTF",
                        "query_end": 293,
                        "query_start": 4,
                        "sbjct": "EFVFIY",
                        "sbjct_end": 342,
                        "sbjct_start": 6
                    }
                ],
                "length": 343,
                "title": "SouthBay"
            }
        ],
        "query": "BHAV"
    }, sort_keys=True)
]) + '\n'


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
            self.assertEqual('DIAMOND', reader.application)
            self.assertEqual(
                {
                    'reference': (
                        'Buchfink et al., Fast and Sensitive Protein '
                        'Alignment using DIAMOND, Nature Methods, 12, '
                        '59-60 (2015)'),
                    'task': 'blastx',
                    'version': 'v0.8.23',
                },
                reader.params)

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

    def testSaveAsJSON(self):
        """
        A DiamondTabularFormatReader must be able to save itself as JSON.
        """
        mockOpener = mockOpen(read_data=DIAMOND_RECORDS)
        with patch.object(builtins, 'open', mockOpener):
            reader = DiamondTabularFormatReader('file.txt')
            fp = StringIO()
            reader.saveAsJSON(fp)
            self.assertEqual(DIAMOND_RECORDS_DUMPED, fp.getvalue())

    def testSaveAsJSONBzip2(self):
        """
        A DiamondTabularFormatReader must be able to save itself as bzip2'd
        JSON.
        """
        mockOpener = mockOpen(read_data=DIAMOND_RECORDS)
        with patch.object(builtins, 'open', mockOpener):
            reader = DiamondTabularFormatReader('file.txt')
            data = BytesIO()
            fp = bz2file.BZ2File(data, 'w')
            reader.saveAsJSON(fp, writeBytes=True)
            fp.close()
            self.assertEqual(
                compress(DIAMOND_RECORDS_DUMPED.encode('UTF-8')),
                data.getvalue())


_JSON_RECORDS = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59-60 '
                      '(2015)'),
        'task': 'blastx',
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 5075,
                'hsps': [
                    {
                        'sbjct_start': 1817,
                        'sbjct_end': 1849,
                        'bits': 165.393,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 2.73597e-40,
                        'query_end': 99,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTTTGTTGTGGTTG'),
                    }
                ]
            }
        ],
        'query': 'id1'
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 5075,
                'hsps': [
                    {
                        'sbjct_start': 4074,
                        'sbjct_end': 4106,
                        'bits': 178.016,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 5063,
                'hsps': [
                    {
                        'sbjct_start': 4062,
                        'sbjct_end': 4094,
                        'bits': 123.915,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 5046,
                'hsps': [
                    {
                        'sbjct_start': 4039,
                        'sbjct_end': 4070,
                        'bits': 87.848,
                        'frame': 1,
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'id2'
    },
    {
        'alignments': [],
        'query': 'id3'
    },
    {
        'alignments': [],
        'query': 'id4'
    }

]

_JSON_RECORDS_ONE_MIDDLE = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59-60 '
                      '(2015)'),
        'task': 'blastx',
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 5075,
                'hsps': [
                    {
                        'sbjct_start': 1817,
                        'sbjct_end': 1849,
                        'bits': 165.393,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 2.73597e-40,
                        'query_end': 99,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTTTGTTGTGGTTG'),
                    }
                ]
            }
        ],
        'query': 'id1'
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 5075,
                'hsps': [
                    {
                        'sbjct_start': 4074,
                        'sbjct_end': 4106,
                        'bits': 178.016,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 5063,
                'hsps': [
                    {
                        'sbjct_start': 4062,
                        'sbjct_end': 4094,
                        'bits': 123.915,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 5046,
                'hsps': [
                    {
                        'sbjct_start': 4039,
                        'sbjct_end': 4070,
                        'bits': 87.848,
                        'frame': 1,
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'id2'
    },
    {
        'alignments': [],
        'query': 'id4'
    }

]

_JSON_RECORDS_ONE_END = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59-60 '
                      '(2015)'),
        'task': 'blastx',
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 5075,
                'hsps': [
                    {
                        'sbjct_start': 1817,
                        'sbjct_end': 1849,
                        'bits': 165.393,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 2.73597e-40,
                        'query_end': 99,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTTTGTTGTGGTTG'),
                    }
                ]
            }
        ],
        'query': 'id1'
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 5075,
                'hsps': [
                    {
                        'sbjct_start': 4074,
                        'sbjct_end': 4106,
                        'bits': 178.016,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 5063,
                'hsps': [
                    {
                        'sbjct_start': 4062,
                        'sbjct_end': 4094,
                        'bits': 123.915,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 5046,
                'hsps': [
                    {
                        'sbjct_start': 4039,
                        'sbjct_end': 4070,
                        'bits': 87.848,
                        'frame': 1,
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'id2'
    },
    {
        'alignments': [],
        'query': 'id3'
    },

]

_JSON_RECORDS_ONE_START = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59-60 '
                      '(2015)'),
        'task': 'blastx',
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 5075,
                'hsps': [
                    {
                        'sbjct_start': 4074,
                        'sbjct_end': 4106,
                        'bits': 178.016,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 5063,
                'hsps': [
                    {
                        'sbjct_start': 4062,
                        'sbjct_end': 4094,
                        'bits': 123.915,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 5046,
                'hsps': [
                    {
                        'sbjct_start': 4039,
                        'sbjct_end': 4070,
                        'bits': 87.848,
                        'frame': 1,
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'id2'
    },
    {
        'alignments': [],
        'query': 'id3'
    },
    {
        'alignments': [],
        'query': 'id4'
    }

]

_JSON_RECORDS_TWO_END = [
    {
        'application': 'DIAMOND',
        'version': 'v0.8.23',
        'reference': ('Buchfink et al., Fast and Sensitive Protein '
                      'Alignment using DIAMOND, Nature Methods, 12, 59-60 '
                      '(2015)'),
        'task': 'blastx',
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 5075,
                'hsps': [
                    {
                        'sbjct_start': 1817,
                        'sbjct_end': 1849,
                        'bits': 165.393,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 2.73597e-40,
                        'query_end': 99,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT'
                                  'TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT'
                                  'GTTTTGTTGTGGTTG'),
                    }
                ]
            }
        ],
        'query': 'id1'
    },
    {
        'alignments': [
            {
                'title': 'gi|9629198|ref|NC_001781.1| Human RSV',
                'length': 5075,
                'hsps': [
                    {
                        'sbjct_start': 4074,
                        'sbjct_end': 4106,
                        'bits': 178.016,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 4.33545e-44,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9629367|ref|NC_001803.1| RSV',
                'length': 5063,
                'hsps': [
                    {
                        'sbjct_start': 4062,
                        'sbjct_end': 4094,
                        'bits': 123.915,
                        'frame': 1,
                        'query_start': 1,
                        'expect': 8.37678e-28,
                        'query_end': 101,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGT'
                                  'CCTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCT'
                                  'ATAATTATACCACTGGC'),
                    }
                ]
            },
            {
                'title': 'gi|9631267|ref|NC_001989.1| Bovine RSV',
                'length': 5046,
                'hsps': [
                    {
                        'sbjct_start': 4039,
                        'sbjct_end': 4070,
                        'bits': 87.848,
                        'frame': 1,
                        'query_start': 2,
                        'expect': 6.03169e-17,
                        'query_end': 98,
                        'sbjct': ('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'),
                        'query': ('TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC'
                                  'CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA'
                                  'TAATTATACCACT'),
                    }
                ]
            }
        ],
        'query': 'id2'
    },

]


def _recordsToStr(records):
    """
    Convert a list of DIAMOND JSON records to a string.

    @param records: A C{list} of C{dict}s as would be found in our per-line
        JSON conversion of DIAMOND's tabular output.
    @return: A C{str} suitable for use in a test simulating reading input
        containing our per-line JSON.
    """
    return '\n'.join(dumps(record) for record in records) + '\n'


JSON = _recordsToStr(_JSON_RECORDS)
JSON_ONE_MIDDLE = _recordsToStr(_JSON_RECORDS_ONE_MIDDLE)
JSON_ONE_END = _recordsToStr(_JSON_RECORDS_ONE_END)
JSON_ONE_START = _recordsToStr(_JSON_RECORDS_ONE_START)
JSON_TWO_END = _recordsToStr(_JSON_RECORDS_TWO_END)


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
                'id1',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'id2',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'id3',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'id4',
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
                'id1',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'id2',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'id3',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'id4',
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
        the last read has no matches.  (That read will not be examined by the
        JSONRecordsReader.)
        """
        reads = Reads([
            AARead(
                'id1',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'id2',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'id3',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'id4',
                'GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT'
                'AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC'),
        ])

        mockOpener = mockOpen(read_data=JSON_ONE_END)
        with patch.object(builtins, 'open', mockOpener):
            reader = JSONRecordsReader('file.json')
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(3, len(alignments))

    def testCorrectNumberOfAlignmentsMatchMissingStart(self):
        """
        A JSONRecordsReader must return the expected number of alignments, if
        the first read has no matches.
        """
        reads = Reads([
            AARead(
                'id1',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'id2',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'id3',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'id4',
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
        two reads at the end don't have any matches. (Those reads will not be
        examined by the JSONRecordsReader.)
        """
        reads = Reads([
            AARead(
                'id1',
                'AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG'
                'GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG'),
            AARead(
                'id2',
                'TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC'
                'ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC'),
            AARead(
                'id3',
                'ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA'
                'GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT'),
            AARead(
                'id4',
                'GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT'
                'AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC'),
        ])

        mockOpener = mockOpen(read_data=JSON_TWO_END)
        with patch.object(builtins, 'open', mockOpener):
            reader = JSONRecordsReader('file.json')
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(2, len(alignments))