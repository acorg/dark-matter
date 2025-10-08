import builtins
from bz2 import compress
from io import BytesIO, StringIO
from json import dumps
from unittest import TestCase
from unittest.mock import mock_open, patch

import bz2file

from dark.diamond.conversion import (
    DiamondTabularFormat,
    DiamondTabularFormatReader,
    JSONRecordsReader,
)
from dark.reads import AARead, Reads

# The 15 fields expected in the DIAMOND output we parse are:
#
# qtitle, stitle, bitscore, evalue, qframe, qseq, qstart, qend, sseq, sstart,
# send, slen, btop, nident, positive
#
# See the --outfmt section of 'diamond help' for detail on these directives.
#
# Note that the fields below must be separated by TABs.
DIAMOND_RECORDS = """\
ACC94	INSV	29.6	0.003	1	EFII	178	295	SSSEV	175	285	295	4	0	1
ACC94	CASV	28.1	0.008	1	KLL	7	37	ITRV	9	39	300	3	1	2
ACC94	Golden	28.1	0.009	1	IKSKL	7	35	EETSR	9	37	293	5	2	3
ACC94	Golden	23.5	0.21	1	TIMSVV	177	240	DDMV	179	235	293	6	3	4
ACC94	HCoV	25.0	0.084	1	LHVNYL	1	203	DEELKA	2	210	290	6	4	5
ACC94	HCoV	18.5	9.1	1	SEIICE	226	257	VETVAQ	20	45	290	9	5	6
ACC94	FERV	24.6	0.11	1	YSCFT	176	276	LGKRMFC	152	243	270	10	6	7
AKAV	AKAV	634	0.0	1	GEPFSV	1	306	NIYGEP	1	306	306	8	7	8
AKAV	WYOV	401	7e-143	1	PFSVYG	1	306	GEPMS	1	294	294	8	8	9
BHAV	TAIV	28.1	0.008	1	PKELHG	14	118	SLKSKE	15	131	307	8	9	10
BHAV	SBAY	28.1	0.009	1	CRPTF	4	293	EFVFIY	6	342	343	5	10	11
"""

# The 13 fields expected in the DIAMOND output before we added identities
# and positives (in version 2.0.3) were:
#
# qtitle, stitle, bitscore, evalue, qframe, qseq, qstart, qend, sseq, sstart,
# send, slen, btop
#
# See the --outfmt section of 'diamond help' for detail on these directives.
#
# Note that the fields below must be separated by TABs.
DIAMOND_RECORDS_WITHOUT_NIDENT_AND_POSITIVE = """\
ACC94	INSV	29.6	0.003	1	EFII	178	295	SSSEV	175	285	295	4
ACC94	CASV	28.1	0.008	1	KLL	7	37	ITRV	9	39	300	3
ACC94	Golden	28.1	0.009	1	IKSKL	7	35	EETSR	9	37	293	5
ACC94	Golden	23.5	0.21	1	TIMSVV	177	240	DDMV	179	235	293	6
ACC94	HCoV	25.0	0.084	1	LHVNYL	1	203	DEELKA	2	210	290	6
ACC94	HCoV	18.5	9.1	1	SEIICE	226	257	VETVAQ	20	45	290	9
ACC94	FERV	24.6	0.11	1	YSCFT	176	276	LGKRMFC	152	243	270	10
AKAV	AKAV	634	0.0	1	GEPFSV	1	306	NIYGEP	1	306	306	8
AKAV	WYOV	401	7e-143	1	PFSVYG	1	306	GEPMS	1	294	294	8
BHAV	TAIV	28.1	0.008	1	PKELHG	14	118	SLKSKE	15	131	307	8
BHAV	SBAY	28.1	0.009	1	CRPTF	4	293	EFVFIY	6	342	343	5
"""

# The 17 fields expected in the DIAMOND output with nident, pident (number
# and percent identical) and positive & ppos (number and percent positive).
#
# qtitle, stitle, bitscore, evalue, qframe, qseq, qstart, qend, sseq, sstart,
# send, slen, btop, nident, pident, positive, ppos
#
# See the --outfmt section of 'diamond help' for detail on these directives.
#
# Note that the fields below must be separated by TABs.
DIAMOND_RECORDS_WITH_PERCENTS = """\
ACC94	INSV	29.6	0.003	1	EFII	178	295	SSSEV	175	285	295	4	2	5.0	1	4.0
ACC94	CASV	28.1	0.008	1	KLL	7	37	ITRV	9	39	300	3	3	10.0	2	9.0
ACC94	Golden	28.1	0.009	1	IKSKL	7	35	EETSR	9	37	293	5	4	15.0	3	14.0
ACC94	Golden	23.5	0.21	1	TIMSVV	177	240	DDMV	179	235	293	6	5	20.0	4	19.0
ACC94	HCoV	25.0	0.084	1	LHVNYL	1	203	DEELKA	2	210	290	6	6	25.0	5	24.0
ACC94	HCoV	18.5	9.1	1	SEIICE	226	257	VETVAQ	20	45	290	9	7	30.0	6	29.0
ACC94	FERV	24.6	0.11	1	YSCFT	176	276	LGKRMFC	152	243	270	10	8	35.0	7	34.0
AKAV	AKAV	634	0.0	1	GEPFSV	1	306	NIYGEP	1	306	306	8	9	40.0	8	39.0
AKAV	WYOV	401	7e-143	1	PFSVYG	1	306	GEPMS	1	294	294	8	10	45.0	9	44.0
BHAV	TAIV	28.1	0.008	1	PKELHG	14	118	SLKSKE	15	131	307	8	11	50.0	10	49.0
BHAV	SBAY	28.1	0.009	1	CRPTF	4	293	EFVFIY	6	342	343	5	12	55.0	11	54.0
"""

DIAMOND_RECORD_WITH_SPACES_IN_TITLES = """\
ACC 94	IN SV	29.6	0.003	1	EFII	178	295	SSSEV	175	285	295	4
"""

DIAMOND_RECORDS_DUMPED = (
    "\n".join(
        [
            dumps(
                {
                    "application": "DIAMOND",
                    "reference": (
                        "Buchfink et al., Fast and Sensitive "
                        "Protein Alignment using DIAMOND, Nature Methods, "
                        "12, 59-60 (2015)"
                    ),
                    "task": "blastx",
                    "version": "v0.8.23",
                },
                sort_keys=True,
            ),
            dumps(
                {
                    "alignments": [
                        {
                            "hsps": [
                                {
                                    "bits": 29.6,
                                    "btop": "4",
                                    "expect": 0.003,
                                    "frame": 1,
                                    "identicalCount": 0,
                                    "percentIdentical": None,
                                    "positiveCount": 1,
                                    "percentPositive": None,
                                    "query": "EFII",
                                    "query_end": 295,
                                    "query_start": 178,
                                    "sbjct": "SSSEV",
                                    "sbjct_end": 285,
                                    "sbjct_start": 175,
                                }
                            ],
                            "length": 295,
                            "title": "INSV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 28.1,
                                    "btop": "3",
                                    "expect": 0.008,
                                    "frame": 1,
                                    "identicalCount": 1,
                                    "percentIdentical": None,
                                    "positiveCount": 2,
                                    "percentPositive": None,
                                    "query": "KLL",
                                    "query_end": 37,
                                    "query_start": 7,
                                    "sbjct": "ITRV",
                                    "sbjct_end": 39,
                                    "sbjct_start": 9,
                                }
                            ],
                            "length": 300,
                            "title": "CASV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 28.1,
                                    "btop": "5",
                                    "expect": 0.009,
                                    "frame": 1,
                                    "identicalCount": 2,
                                    "percentIdentical": None,
                                    "positiveCount": 3,
                                    "percentPositive": None,
                                    "query": "IKSKL",
                                    "query_end": 35,
                                    "query_start": 7,
                                    "sbjct": "EETSR",
                                    "sbjct_end": 37,
                                    "sbjct_start": 9,
                                },
                                {
                                    "bits": 23.5,
                                    "btop": "6",
                                    "expect": 0.21,
                                    "frame": 1,
                                    "identicalCount": 3,
                                    "percentIdentical": None,
                                    "positiveCount": 4,
                                    "percentPositive": None,
                                    "query": "TIMSVV",
                                    "query_end": 240,
                                    "query_start": 177,
                                    "sbjct": "DDMV",
                                    "sbjct_end": 235,
                                    "sbjct_start": 179,
                                },
                            ],
                            "length": 293,
                            "title": "Golden",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 25.0,
                                    "btop": "6",
                                    "expect": 0.084,
                                    "frame": 1,
                                    "identicalCount": 4,
                                    "percentIdentical": None,
                                    "positiveCount": 5,
                                    "percentPositive": None,
                                    "query": "LHVNYL",
                                    "query_end": 203,
                                    "query_start": 1,
                                    "sbjct": "DEELKA",
                                    "sbjct_end": 210,
                                    "sbjct_start": 2,
                                },
                                {
                                    "bits": 18.5,
                                    "btop": "9",
                                    "expect": 9.1,
                                    "frame": 1,
                                    "identicalCount": 5,
                                    "percentIdentical": None,
                                    "positiveCount": 6,
                                    "percentPositive": None,
                                    "query": "SEIICE",
                                    "query_end": 257,
                                    "query_start": 226,
                                    "sbjct": "VETVAQ",
                                    "sbjct_end": 45,
                                    "sbjct_start": 20,
                                },
                            ],
                            "length": 290,
                            "title": "HCoV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 24.6,
                                    "btop": "10",
                                    "expect": 0.11,
                                    "frame": 1,
                                    "identicalCount": 6,
                                    "percentIdentical": None,
                                    "positiveCount": 7,
                                    "percentPositive": None,
                                    "query": "YSCFT",
                                    "query_end": 276,
                                    "query_start": 176,
                                    "sbjct": "LGKRMFC",
                                    "sbjct_end": 243,
                                    "sbjct_start": 152,
                                }
                            ],
                            "length": 270,
                            "title": "FERV",
                        },
                    ],
                    "query": "ACC94",
                },
                sort_keys=True,
            ),
            dumps(
                {
                    "alignments": [
                        {
                            "hsps": [
                                {
                                    "bits": 634.0,
                                    "btop": "8",
                                    "expect": 0.0,
                                    "frame": 1,
                                    "identicalCount": 7,
                                    "percentIdentical": None,
                                    "positiveCount": 8,
                                    "percentPositive": None,
                                    "query": "GEPFSV",
                                    "query_end": 306,
                                    "query_start": 1,
                                    "sbjct": "NIYGEP",
                                    "sbjct_end": 306,
                                    "sbjct_start": 1,
                                }
                            ],
                            "length": 306,
                            "title": "AKAV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 401.0,
                                    "btop": "8",
                                    "expect": 7e-143,
                                    "frame": 1,
                                    "identicalCount": 8,
                                    "percentIdentical": None,
                                    "positiveCount": 9,
                                    "percentPositive": None,
                                    "query": "PFSVYG",
                                    "query_end": 306,
                                    "query_start": 1,
                                    "sbjct": "GEPMS",
                                    "sbjct_end": 294,
                                    "sbjct_start": 1,
                                }
                            ],
                            "length": 294,
                            "title": "WYOV",
                        },
                    ],
                    "query": "AKAV",
                },
                sort_keys=True,
            ),
            dumps(
                {
                    "alignments": [
                        {
                            "hsps": [
                                {
                                    "bits": 28.1,
                                    "btop": "8",
                                    "expect": 0.008,
                                    "frame": 1,
                                    "identicalCount": 9,
                                    "percentIdentical": None,
                                    "positiveCount": 10,
                                    "percentPositive": None,
                                    "query": "PKELHG",
                                    "query_end": 118,
                                    "query_start": 14,
                                    "sbjct": "SLKSKE",
                                    "sbjct_end": 131,
                                    "sbjct_start": 15,
                                }
                            ],
                            "length": 307,
                            "title": "TAIV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 28.1,
                                    "btop": "5",
                                    "expect": 0.009,
                                    "frame": 1,
                                    "identicalCount": 10,
                                    "percentIdentical": None,
                                    "positiveCount": 11,
                                    "percentPositive": None,
                                    "query": "CRPTF",
                                    "query_end": 293,
                                    "query_start": 4,
                                    "sbjct": "EFVFIY",
                                    "sbjct_end": 342,
                                    "sbjct_start": 6,
                                }
                            ],
                            "length": 343,
                            "title": "SBAY",
                        },
                    ],
                    "query": "BHAV",
                },
                sort_keys=True,
            ),
        ]
    )
    + "\n"
)


DIAMOND_RECORDS_WITH_PERCENTS_DUMPED = (
    "\n".join(
        [
            dumps(
                {
                    "application": "DIAMOND",
                    "reference": (
                        "Buchfink et al., Fast and Sensitive "
                        "Protein Alignment using DIAMOND, Nature Methods, "
                        "12, 59-60 (2015)"
                    ),
                    "task": "blastx",
                    "version": "v0.8.23",
                },
                sort_keys=True,
            ),
            dumps(
                {
                    "alignments": [
                        {
                            "hsps": [
                                {
                                    "bits": 29.6,
                                    "btop": "4",
                                    "expect": 0.003,
                                    "frame": 1,
                                    "identicalCount": 2,
                                    "percentIdentical": 5.0,
                                    "percentPositive": 4.0,
                                    "positiveCount": 1,
                                    "query": "EFII",
                                    "query_end": 295,
                                    "query_start": 178,
                                    "sbjct": "SSSEV",
                                    "sbjct_end": 285,
                                    "sbjct_start": 175,
                                }
                            ],
                            "length": 295,
                            "title": "INSV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 28.1,
                                    "btop": "3",
                                    "expect": 0.008,
                                    "frame": 1,
                                    "identicalCount": 3,
                                    "percentIdentical": 10.0,
                                    "percentPositive": 9.0,
                                    "positiveCount": 2,
                                    "query": "KLL",
                                    "query_end": 37,
                                    "query_start": 7,
                                    "sbjct": "ITRV",
                                    "sbjct_end": 39,
                                    "sbjct_start": 9,
                                }
                            ],
                            "length": 300,
                            "title": "CASV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 28.1,
                                    "btop": "5",
                                    "expect": 0.009,
                                    "frame": 1,
                                    "identicalCount": 4,
                                    "percentIdentical": 15.0,
                                    "percentPositive": 14.0,
                                    "positiveCount": 3,
                                    "query": "IKSKL",
                                    "query_end": 35,
                                    "query_start": 7,
                                    "sbjct": "EETSR",
                                    "sbjct_end": 37,
                                    "sbjct_start": 9,
                                },
                                {
                                    "bits": 23.5,
                                    "btop": "6",
                                    "expect": 0.21,
                                    "frame": 1,
                                    "identicalCount": 5,
                                    "percentIdentical": 20.0,
                                    "percentPositive": 19.0,
                                    "positiveCount": 4,
                                    "query": "TIMSVV",
                                    "query_end": 240,
                                    "query_start": 177,
                                    "sbjct": "DDMV",
                                    "sbjct_end": 235,
                                    "sbjct_start": 179,
                                },
                            ],
                            "length": 293,
                            "title": "Golden",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 25.0,
                                    "btop": "6",
                                    "expect": 0.084,
                                    "frame": 1,
                                    "identicalCount": 6,
                                    "percentIdentical": 25.0,
                                    "percentPositive": 24.0,
                                    "positiveCount": 5,
                                    "query": "LHVNYL",
                                    "query_end": 203,
                                    "query_start": 1,
                                    "sbjct": "DEELKA",
                                    "sbjct_end": 210,
                                    "sbjct_start": 2,
                                },
                                {
                                    "bits": 18.5,
                                    "btop": "9",
                                    "expect": 9.1,
                                    "frame": 1,
                                    "identicalCount": 7,
                                    "percentIdentical": 30.0,
                                    "percentPositive": 29.0,
                                    "positiveCount": 6,
                                    "query": "SEIICE",
                                    "query_end": 257,
                                    "query_start": 226,
                                    "sbjct": "VETVAQ",
                                    "sbjct_end": 45,
                                    "sbjct_start": 20,
                                },
                            ],
                            "length": 290,
                            "title": "HCoV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 24.6,
                                    "btop": "10",
                                    "expect": 0.11,
                                    "frame": 1,
                                    "identicalCount": 8,
                                    "percentIdentical": 35.0,
                                    "percentPositive": 34.0,
                                    "positiveCount": 7,
                                    "query": "YSCFT",
                                    "query_end": 276,
                                    "query_start": 176,
                                    "sbjct": "LGKRMFC",
                                    "sbjct_end": 243,
                                    "sbjct_start": 152,
                                }
                            ],
                            "length": 270,
                            "title": "FERV",
                        },
                    ],
                    "query": "ACC94",
                },
                sort_keys=True,
            ),
            dumps(
                {
                    "alignments": [
                        {
                            "hsps": [
                                {
                                    "bits": 634.0,
                                    "btop": "8",
                                    "expect": 0.0,
                                    "frame": 1,
                                    "identicalCount": 9,
                                    "percentIdentical": 40.0,
                                    "percentPositive": 39.0,
                                    "positiveCount": 8,
                                    "query": "GEPFSV",
                                    "query_end": 306,
                                    "query_start": 1,
                                    "sbjct": "NIYGEP",
                                    "sbjct_end": 306,
                                    "sbjct_start": 1,
                                }
                            ],
                            "length": 306,
                            "title": "AKAV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 401.0,
                                    "btop": "8",
                                    "expect": 7e-143,
                                    "frame": 1,
                                    "identicalCount": 10,
                                    "percentIdentical": 45.0,
                                    "percentPositive": 44.0,
                                    "positiveCount": 9,
                                    "query": "PFSVYG",
                                    "query_end": 306,
                                    "query_start": 1,
                                    "sbjct": "GEPMS",
                                    "sbjct_end": 294,
                                    "sbjct_start": 1,
                                }
                            ],
                            "length": 294,
                            "title": "WYOV",
                        },
                    ],
                    "query": "AKAV",
                },
                sort_keys=True,
            ),
            dumps(
                {
                    "alignments": [
                        {
                            "hsps": [
                                {
                                    "bits": 28.1,
                                    "btop": "8",
                                    "expect": 0.008,
                                    "frame": 1,
                                    "identicalCount": 11,
                                    "percentIdentical": 50.0,
                                    "percentPositive": 49.0,
                                    "positiveCount": 10,
                                    "query": "PKELHG",
                                    "query_end": 118,
                                    "query_start": 14,
                                    "sbjct": "SLKSKE",
                                    "sbjct_end": 131,
                                    "sbjct_start": 15,
                                }
                            ],
                            "length": 307,
                            "title": "TAIV",
                        },
                        {
                            "hsps": [
                                {
                                    "bits": 28.1,
                                    "btop": "5",
                                    "expect": 0.009,
                                    "frame": 1,
                                    "identicalCount": 12,
                                    "percentIdentical": 55.0,
                                    "percentPositive": 54.0,
                                    "positiveCount": 11,
                                    "query": "CRPTF",
                                    "query_end": 293,
                                    "query_start": 4,
                                    "sbjct": "EFVFIY",
                                    "sbjct_end": 342,
                                    "sbjct_start": 6,
                                }
                            ],
                            "length": 343,
                            "title": "SBAY",
                        },
                    ],
                    "query": "BHAV",
                },
                sort_keys=True,
            ),
        ]
    )
    + "\n"
)


class TestDiamondTabularFormatReader(TestCase):
    """
    Test the DiamondTabularFormatReader class.
    """

    def testDiamondParams(self):
        """
        When a DIAMOND file has been read, its parameters must be present
        in the reader instance.
        """
        mockOpener = mock_open(read_data=DIAMOND_RECORDS)
        with patch.object(builtins, "open", mockOpener):
            reader = DiamondTabularFormatReader("file.txt")
            list(reader.records())
            self.assertEqual("DIAMOND", reader.application)
            self.assertEqual(
                {
                    "application": "DIAMOND",
                    "reference": (
                        "Buchfink et al., Fast and Sensitive Protein "
                        "Alignment using DIAMOND, Nature Methods, 12, "
                        "59-60 (2015)"
                    ),
                    "task": "blastx",
                    "version": "v0.8.23",
                },
                reader.params,
            )

    def testDiamondInput(self):
        """
        Test conversion of a chunk of DIAMOND output.
        """
        mockOpener = mock_open(read_data=DIAMOND_RECORDS)
        with patch.object(builtins, "open", mockOpener):
            reader = DiamondTabularFormatReader("file.txt")
            acc94, akav, bhav = list(reader.records())
            self.assertEqual(5, len(acc94["alignments"]))
            self.assertEqual(2, len(akav["alignments"]))
            self.assertEqual(2, len(bhav["alignments"]))

    def testDiamondInputWithoutNidentOrPositives(self):
        """
        Test conversion of a chunk of DIAMOND output that does not contain the
        nident, positives, or pident fields.
        """
        mockOpener = mock_open(read_data=DIAMOND_RECORDS_WITHOUT_NIDENT_AND_POSITIVE)
        with patch.object(builtins, "open", mockOpener):
            reader = DiamondTabularFormatReader("file.txt")
            acc94, akav, bhav = list(reader.records())
            for record in acc94, akav, bhav:
                for alignment in record["alignments"]:
                    for hsp in alignment["hsps"]:
                        self.assertIs(None, hsp["identicalCount"])
                        self.assertIs(None, hsp["positiveCount"])
                        self.assertIs(None, hsp["percentIdentical"])

    def testSaveAsJSON(self):
        """
        A DiamondTabularFormatReader must be able to save itself as JSON.
        """
        mockOpener = mock_open(read_data=DIAMOND_RECORDS)
        with patch.object(builtins, "open", mockOpener):
            reader = DiamondTabularFormatReader("file.txt")
            fp = StringIO()
            reader.saveAsJSON(fp)
            self.maxDiff = None
            self.assertEqual(DIAMOND_RECORDS_DUMPED, fp.getvalue())

    def testSaveAsJSONWithPercentIdentical(self):
        """
        A DiamondTabularFormatReader must be able to save itself as JSON
        when the percentIdentical field is present.
        """
        mockOpener = mock_open(read_data=DIAMOND_RECORDS_WITH_PERCENTS)
        with patch.object(builtins, "open", mockOpener):
            reader = DiamondTabularFormatReader("file.txt")
            fp = StringIO()
            reader.saveAsJSON(fp)
            self.maxDiff = None
            self.assertEqual(DIAMOND_RECORDS_WITH_PERCENTS_DUMPED, fp.getvalue())

    def testSaveAsJSONBzip2(self):
        """
        A DiamondTabularFormatReader must be able to save itself as bzip2'd
        JSON.
        """
        mockOpener = mock_open(read_data=DIAMOND_RECORDS)
        with patch.object(builtins, "open", mockOpener):
            reader = DiamondTabularFormatReader("file.txt")
            data = BytesIO()
            fp = bz2file.BZ2File(data, "w")
            reader.saveAsJSON(fp, writeBytes=True)
            fp.close()
            self.assertEqual(
                compress(DIAMOND_RECORDS_DUMPED.encode("UTF-8")), data.getvalue()
            )

    def testSaveAsJSONBzip2WithPercentIdentical(self):
        """
        A DiamondTabularFormatReader must be able to save itself as bzip2'd
        JSON when the percentIdentical field is present.
        """
        mockOpener = mock_open(read_data=DIAMOND_RECORDS_WITH_PERCENTS)
        with patch.object(builtins, "open", mockOpener):
            reader = DiamondTabularFormatReader("file.txt")
            data = BytesIO()
            fp = bz2file.BZ2File(data, "w")
            reader.saveAsJSON(fp, writeBytes=True)
            fp.close()
            self.assertEqual(
                compress(DIAMOND_RECORDS_WITH_PERCENTS_DUMPED.encode("UTF-8")),
                data.getvalue(),
            )

    def testSpacesMustBePreserved(self):
        """
        If there are spaces in the query title or subject titles, the spaces
        must be preserved.
        """
        mockOpener = mock_open(read_data=DIAMOND_RECORD_WITH_SPACES_IN_TITLES)
        with patch.object(builtins, "open", mockOpener):
            reader = DiamondTabularFormatReader("file.txt")
            acc94 = list(reader.records())
            self.assertEqual("ACC 94", acc94[0]["query"])
            self.assertEqual("IN SV", acc94[0]["alignments"][0]["title"])


_JSON_RECORDS = [
    {
        "application": "DIAMOND",
        "version": "v0.8.23",
        "reference": (
            "Buchfink et al., Fast and Sensitive Protein "
            "Alignment using DIAMOND, Nature Methods, 12, 59-60 "
            "(2015)"
        ),
        "task": "blastx",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 5075,
                "hsps": [
                    {
                        "sbjct_start": 1817,
                        "sbjct_end": 1849,
                        "bits": 165.393,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 2.73597e-40,
                        "query_end": 99,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                        "query": (
                            "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT"
                            "TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT"
                            "GTTTTGTTGTGGTTG"
                        ),
                    }
                ],
            }
        ],
        "query": "id1",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 5075,
                "hsps": [
                    {
                        "sbjct_start": 4074,
                        "sbjct_end": 4106,
                        "bits": 178.016,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 4.33545e-44,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5063,
                "hsps": [
                    {
                        "sbjct_start": 4062,
                        "sbjct_end": 4094,
                        "bits": 123.915,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 8.37678e-28,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5046,
                "hsps": [
                    {
                        "sbjct_start": 4039,
                        "sbjct_end": 4070,
                        "bits": 87.848,
                        "btop": "",
                        "frame": 1,
                        "query_start": 2,
                        "expect": 6.03169e-17,
                        "query_end": 98,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                        "query": (
                            "TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC"
                            "CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA"
                            "TAATTATACCACT"
                        ),
                    }
                ],
            },
        ],
        "query": "id2",
    },
    {"alignments": [], "query": "id3"},
    {"alignments": [], "query": "id4"},
]

_JSON_RECORDS_ONE_MIDDLE = [
    {
        "application": "DIAMOND",
        "version": "v0.8.23",
        "reference": (
            "Buchfink et al., Fast and Sensitive Protein "
            "Alignment using DIAMOND, Nature Methods, 12, 59-60 "
            "(2015)"
        ),
        "task": "blastx",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 5075,
                "hsps": [
                    {
                        "sbjct_start": 1817,
                        "sbjct_end": 1849,
                        "bits": 165.393,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 2.73597e-40,
                        "query_end": 99,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                        "query": (
                            "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT"
                            "TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT"
                            "GTTTTGTTGTGGTTG"
                        ),
                    }
                ],
            }
        ],
        "query": "id1",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 5075,
                "hsps": [
                    {
                        "sbjct_start": 4074,
                        "sbjct_end": 4106,
                        "bits": 178.016,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 4.33545e-44,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5063,
                "hsps": [
                    {
                        "sbjct_start": 4062,
                        "sbjct_end": 4094,
                        "bits": 123.915,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 8.37678e-28,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5046,
                "hsps": [
                    {
                        "sbjct_start": 4039,
                        "sbjct_end": 4070,
                        "bits": 87.848,
                        "btop": "",
                        "frame": 1,
                        "query_start": 2,
                        "expect": 6.03169e-17,
                        "query_end": 98,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                        "query": (
                            "TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC"
                            "CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA"
                            "TAATTATACCACT"
                        ),
                    }
                ],
            },
        ],
        "query": "id2",
    },
    {"alignments": [], "query": "id4"},
]

_JSON_RECORDS_ONE_END = [
    {
        "application": "DIAMOND",
        "version": "v0.8.23",
        "reference": (
            "Buchfink et al., Fast and Sensitive Protein "
            "Alignment using DIAMOND, Nature Methods, 12, 59-60 "
            "(2015)"
        ),
        "task": "blastx",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 5075,
                "hsps": [
                    {
                        "sbjct_start": 1817,
                        "sbjct_end": 1849,
                        "bits": 165.393,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 2.73597e-40,
                        "query_end": 99,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                        "query": (
                            "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT"
                            "TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT"
                            "GTTTTGTTGTGGTTG"
                        ),
                    }
                ],
            }
        ],
        "query": "id1",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 5075,
                "hsps": [
                    {
                        "sbjct_start": 4074,
                        "sbjct_end": 4106,
                        "bits": 178.016,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 4.33545e-44,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5063,
                "hsps": [
                    {
                        "sbjct_start": 4062,
                        "sbjct_end": 4094,
                        "bits": 123.915,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 8.37678e-28,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5046,
                "hsps": [
                    {
                        "sbjct_start": 4039,
                        "sbjct_end": 4070,
                        "bits": 87.848,
                        "btop": "",
                        "frame": 1,
                        "query_start": 2,
                        "expect": 6.03169e-17,
                        "query_end": 98,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                        "query": (
                            "TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC"
                            "CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA"
                            "TAATTATACCACT"
                        ),
                    }
                ],
            },
        ],
        "query": "id2",
    },
    {"alignments": [], "query": "id3"},
]

_JSON_RECORDS_ONE_START = [
    {
        "application": "DIAMOND",
        "version": "v0.8.23",
        "reference": (
            "Buchfink et al., Fast and Sensitive Protein "
            "Alignment using DIAMOND, Nature Methods, 12, 59-60 "
            "(2015)"
        ),
        "task": "blastx",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 5075,
                "hsps": [
                    {
                        "sbjct_start": 4074,
                        "sbjct_end": 4106,
                        "bits": 178.016,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 4.33545e-44,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5063,
                "hsps": [
                    {
                        "sbjct_start": 4062,
                        "sbjct_end": 4094,
                        "bits": 123.915,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 8.37678e-28,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5046,
                "hsps": [
                    {
                        "sbjct_start": 4039,
                        "sbjct_end": 4070,
                        "bits": 87.848,
                        "btop": "",
                        "frame": 1,
                        "query_start": 2,
                        "expect": 6.03169e-17,
                        "query_end": 98,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                        "query": (
                            "TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC"
                            "CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA"
                            "TAATTATACCACT"
                        ),
                    }
                ],
            },
        ],
        "query": "id2",
    },
    {"alignments": [], "query": "id3"},
    {"alignments": [], "query": "id4"},
]

_JSON_RECORDS_TWO_END = [
    {
        "application": "DIAMOND",
        "version": "v0.8.23",
        "reference": (
            "Buchfink et al., Fast and Sensitive Protein "
            "Alignment using DIAMOND, Nature Methods, 12, 59-60 "
            "(2015)"
        ),
        "task": "blastx",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 5075,
                "hsps": [
                    {
                        "sbjct_start": 1817,
                        "sbjct_end": 1849,
                        "bits": 165.393,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 2.73597e-40,
                        "query_end": 99,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                        "query": (
                            "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGT"
                            "TTTCGGGGGTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGT"
                            "GTTTTGTTGTGGTTG"
                        ),
                    }
                ],
            }
        ],
        "query": "id1 1",
    },
    {
        "alignments": [
            {
                "title": "gi|9629198|ref|NC_001781.1| Human RSV",
                "length": 5075,
                "hsps": [
                    {
                        "sbjct_start": 4074,
                        "sbjct_end": 4106,
                        "bits": 178.016,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 4.33545e-44,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5063,
                "hsps": [
                    {
                        "sbjct_start": 4062,
                        "sbjct_end": 4094,
                        "bits": 123.915,
                        "btop": "",
                        "frame": 1,
                        "query_start": 1,
                        "expect": 8.37678e-28,
                        "query_end": 101,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
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
                "length": 5046,
                "hsps": [
                    {
                        "sbjct_start": 4039,
                        "sbjct_end": 4070,
                        "bits": 87.848,
                        "btop": "",
                        "frame": 1,
                        "query_start": 2,
                        "expect": 6.03169e-17,
                        "query_end": 98,
                        "sbjct": ("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
                        "query": (
                            "TTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTC"
                            "CTCTTTCACCACGAGTTAAACCATTAACATTATATTTTTCTA"
                            "TAATTATACCACT"
                        ),
                    }
                ],
            },
        ],
        "query": "id2 2",
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
    return "\n".join(dumps(record) for record in records) + "\n"


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
        reads = Reads(
            [
                AARead(
                    "id1",
                    "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG"
                    "GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG",
                ),
                AARead(
                    "id2",
                    "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC"
                    "ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC",
                ),
                AARead(
                    "id3",
                    "ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA"
                    "GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT",
                ),
                AARead(
                    "id4",
                    "GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT"
                    "AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC",
                ),
            ]
        )

        mockOpener = mock_open(read_data=JSON)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(alignments))

    def testCorrectNumberOfAlignmentsMatchMissingMiddle(self):
        """
        A JSONRecordsReader must return the expected number of alignments, if
        a match is missing in the middle of the JSON file.
        """
        reads = Reads(
            [
                AARead(
                    "id1",
                    "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG"
                    "GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG",
                ),
                AARead(
                    "id2",
                    "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC"
                    "ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC",
                ),
                AARead(
                    "id3",
                    "ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA"
                    "GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT",
                ),
                AARead(
                    "id4",
                    "GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT"
                    "AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC",
                ),
            ]
        )

        mockOpener = mock_open(read_data=JSON_ONE_MIDDLE)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(alignments))

    def testCorrectNumberOfAlignmentsMatchMissingEnd(self):
        """
        A JSONRecordsReader must return the expected number of alignments, if
        the last read has no matches.  (That read will not be examined by the
        JSONRecordsReader.)
        """
        reads = Reads(
            [
                AARead(
                    "id1",
                    "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG"
                    "GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG",
                ),
                AARead(
                    "id2",
                    "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC"
                    "ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC",
                ),
                AARead(
                    "id3",
                    "ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA"
                    "GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT",
                ),
                AARead(
                    "id4",
                    "GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT"
                    "AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC",
                ),
            ]
        )

        mockOpener = mock_open(read_data=JSON_ONE_END)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(3, len(alignments))

    def testCorrectNumberOfAlignmentsMatchMissingStart(self):
        """
        A JSONRecordsReader must return the expected number of alignments, if
        the first read has no matches.
        """
        reads = Reads(
            [
                AARead(
                    "id1",
                    "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG"
                    "GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG",
                ),
                AARead(
                    "id2",
                    "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC"
                    "ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC",
                ),
                AARead(
                    "id3",
                    "ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA"
                    "GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT",
                ),
                AARead(
                    "id4",
                    "GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT"
                    "AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC",
                ),
            ]
        )

        mockOpener = mock_open(read_data=JSON_ONE_START)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(4, len(alignments))

    def testCorrectNumberOfAlignmentsTwoMatchesMissingEnd(self):
        """
        A JSONRecordsReader must return the expected number of alignments, if
        two reads at the end don't have any matches. (Those reads will not be
        examined by the JSONRecordsReader.)
        """
        reads = Reads(
            [
                AARead(
                    "id1 1",
                    "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG"
                    "GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG",
                ),
                AARead(
                    "id2 2",
                    "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC"
                    "ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC",
                ),
                AARead(
                    "id3 3",
                    "ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA"
                    "GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT",
                ),
                AARead(
                    "id4 4",
                    "GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT"
                    "AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC",
                ),
            ]
        )

        mockOpener = mock_open(read_data=JSON_TWO_END)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            alignments = list(reader.readAlignments(reads))
            self.assertEqual(2, len(alignments))

    def testSpacesMustBePreserved(self):
        """
        A JSONRecordsReader must return the right query and subject titles,
        even if they have spaces.
        """
        reads = Reads(
            [
                AARead(
                    "id1 1",
                    "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG"
                    "GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG",
                ),
                AARead(
                    "id2 2",
                    "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC"
                    "ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC",
                ),
                AARead(
                    "id3 3",
                    "ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA"
                    "GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT",
                ),
                AARead(
                    "id4 4",
                    "GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT"
                    "AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC",
                ),
            ]
        )

        mockOpener = mock_open(read_data=JSON_TWO_END)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            alignment = list(reader.readAlignments(reads))[0]
            self.assertEqual("id1 1", alignment.read.id)

    def testSpaceInReadIdNotInJSONRecord(self):
        """
        A JSONRecordsReader must return the right query and subject titles,
        when the read ids have spaces in them but the titles in the JSON have
        been truncated at the first space (as in the SAM format output of the
        BWA 'mem' command).
        """
        reads = Reads(
            [
                AARead(
                    "id1 1",
                    "AGGGCTCGGATGCTGTGGGTGTTTGTGTGGAGTTGGGTGTGTTTTCGGGG"
                    "GTGGTTGAGTGGAGGGATTGCTGTTGGATTGTGTGTTTTGTTGTGGTTGCG",
                ),
                AARead(
                    "id2 2",
                    "TTTTTCTCCTGCGTAGATGAACCTACCCATGGCTTAGTAGGTCCTCTTTC"
                    "ACCACGAGTTAAACCATTAACATTATATTTTTCTATAATTATACCACTGGC",
                ),
                AARead(
                    "id3 3",
                    "ACCTCCGCCTCCCAGGTTCAAGCAATTCTCCTGCCTTAGCCTCCTGAATA"
                    "GCTGGGATTACAGGTATGCAGGAGGCTAAGGCAGGAGAATTGCTTGAACCT",
                ),
                AARead(
                    "id4 4",
                    "GAGGGTGGAGGTAACTGAGGAAGCAAAGGCTTGGAGACAGGGCCCCTCAT"
                    "AGCCAGTGAGTGCGCCATTTTCTTTGGAGCAATTGGGTGGGGAGATGGGGC",
                ),
            ]
        )

        mockOpener = mock_open(read_data=JSON)
        with patch.object(builtins, "open", mockOpener):
            reader = JSONRecordsReader("file.json")
            alignment = list(reader.readAlignments(reads))[0]
            self.assertEqual("id1 1", alignment.read.id)


class TestDiamondTabularFormatToDicts(TestCase):
    """
    Tests for the diamondTabularFormatToDicts function.
    """

    def testDuplicatesInFieldNameList(self):
        """
        If a field name list that contains duplicates is passed, the
        DiamondTabularFormat __init__ function must raise a ValueError.
        """
        error = "^field names contains duplicated names: a, b\\.$"
        self.assertRaisesRegex(
            ValueError, error, DiamondTabularFormat, ["a", "b", "a", "c", "b"]
        )

    def testTooFewFields(self):
        """
        If an input line does not have enough fields, a ValueError must be
        raised.
        """
        dtf = DiamondTabularFormat(["a", "b", "c"])
        data = StringIO("a\tb\n")
        error = (
            r"^DIAMOND output line had 2 field values \(expected 3\)\. "
            r"The offending input line was 'a\\tb\\n'\."
        )
        self.assertRaisesRegex(
            ValueError, error, list, dtf.diamondTabularFormatToDicts(data)
        )

    def testTooManyFields(self):
        """
        If an input line has too many fields, a ValueError must be raised.
        """
        dtf = DiamondTabularFormat(["a", "b"])
        data = StringIO("a\tb\tc\n")
        error = (
            r"^DIAMOND output line had 3 field values \(expected 2\)\. "
            r"The offending input line was 'a\\tb\\tc\\n'\."
        )
        self.assertRaisesRegex(
            ValueError, error, list, dtf.diamondTabularFormatToDicts(data)
        )

    def testUnknownField(self):
        """
        An unknown field name must result in a returned field name and value
        that are identical to those in the function call and its input string.
        """
        dtf = DiamondTabularFormat(["__blah__"])
        data = StringIO("3.5\n")
        (result,) = list(dtf.diamondTabularFormatToDicts(data))
        self.assertEqual({"__blah__": "3.5"}, result)

    def testConversions(self):
        """
        The fields in input lines must be recognized and converted to their
        correct types.
        """
        fields = [
            "bitscore",
            "evalue",
            "frame",
            "identicalCount",
            "positiveCount",
            "qstart",
            "qend",
            "sstart",
            "send",
            "qseq",
        ]
        data = StringIO(
            ("3.5 1.7 1 7 4 10 12 1 2 ACGT\n3.6 1.8 2 8 5 11 13 2 3 TGCA").replace(
                " ", "\t"
            )
            + "\n"
        )
        dtf = DiamondTabularFormat(fields)
        (result1, result2) = list(dtf.diamondTabularFormatToDicts(data))

        self.assertEqual(
            {
                "bitscore": 3.5,
                "evalue": 1.7,
                "frame": 1,
                "identicalCount": 7,
                "positiveCount": 4,
                "qstart": 10,
                "qend": 12,
                "sstart": 1,
                "send": 2,
                "qseq": "ACGT",
            },
            result1,
        )

        self.assertEqual(
            {
                "bitscore": 3.6,
                "evalue": 1.8,
                "frame": 2,
                "identicalCount": 8,
                "positiveCount": 5,
                "qstart": 11,
                "qend": 13,
                "sstart": 2,
                "send": 3,
                "qseq": "TGCA",
            },
            result2,
        )
