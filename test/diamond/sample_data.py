# Sample DIAMOND parameters.
PARAMS = {
    "application": "DIAMOND",
    "reference": (
        "Buchfink et al., Fast and Sensitive Protein Alignment "
        "using DIAMOND, Nature Methods, 12, 59-60 (2015)"
    ),
    "task": "blastx",
    "version": "v0.8.23",
}

RECORD0 = {
    "query": "id0",
    "alignments": [
        {
            "length": 37000,
            "hsps": [
                {
                    "btop": "",
                    "bits": 20,
                    "sbjct_end": 15400,
                    "expect": 1e-11,
                    "sbjct": "AAAAAAAAAAA",
                    "sbjct_start": 15390,
                    "query": "TACCCTGCGGCCCGCTACGGCTGGTCTCCATTT",
                    "frame": 1,
                    "query_end": 68,
                    "query_start": 36,
                }
            ],
            "title": "gi|887699|gb|DQ37780 Squirrelpox virus 1296/99",
        },
        {
            "length": 38000,
            "hsps": [
                {
                    "btop": "",
                    "bits": 25,
                    "sbjct_end": 12400,
                    "expect": 1e-10,
                    "sbjct": "AAAAAAAAAAA",
                    "sbjct_start": 12390,
                    "query": "TACCCTGCGGCCCGCTACGGCTGGTCTCCATTT",
                    "frame": 1,
                    "query_end": 68,
                    "query_start": 36,
                }
            ],
            "title": "gi|887699|gb|DQ37780 Squirrelpox virus 55",
        },
    ],
}

RECORD1 = {
    "query": "id1",
    "alignments": [
        {
            "length": 35000,
            "hsps": [
                {
                    "btop": "",
                    "bits": 20,
                    "sbjct_end": 11400,
                    "expect": 1e-8,
                    "sbjct": "AAAAAAAAAAA",
                    "sbjct_start": 11390,
                    "query": "TACCCTGCGGCCCGCTACGGCTGGTCTCCATTT",
                    "frame": 1,
                    "query_end": 68,
                    "query_start": 36,
                }
            ],
            "title": "gi|887699|gb|DQ37780 Monkeypox virus 456",
        },
        {
            "length": 35000,
            "hsps": [
                {
                    "btop": "",
                    "bits": 20,
                    "sbjct_end": 10400,
                    "expect": 1e-7,
                    "sbjct": "AAAAA-AAAAAA",
                    "sbjct_start": 10390,
                    "query": "TACCCTGCGGCCCGCTACGGCTGGTCTCCATTT",
                    "frame": 1,
                    "query_end": 68,
                    "query_start": 36,
                }
            ],
            "title": "gi|887699|gb|DQ37780 Mummypox virus 3000 B.C.",
        },
    ],
}

RECORD2 = {
    "query": "id2",
    "alignments": [
        {
            "length": 30000,
            "hsps": [
                {
                    "btop": "",
                    "bits": 20,
                    "sbjct_end": 1400,
                    "expect": 1e-6,
                    "sbjct": "AAAAAAAAAAA",
                    "sbjct_start": 1390,
                    "query": "TACCCTGCGGCCCGCTACGGCTGGTCTCCATTT",
                    "frame": 1,
                    "query_end": 68,
                    "query_start": 36,
                }
            ],
            "title": "gi|887699|gb|DQ37780 Cowpox virus 15",
        }
    ],
}

# Identical to RECORD2, apart from e-value.
RECORD3 = {
    "query": "id3",
    "alignments": [
        {
            "length": 30000,
            "hsps": [
                {
                    "btop": "",
                    "bits": 20,
                    "sbjct_end": 1400,
                    "expect": 1e-5,
                    "sbjct": "AAAAAAAAAAA",
                    "sbjct_start": 1390,
                    "query": "TACCCTGCGGCCCGCTACGGCTGGTCTCCATTT",
                    "frame": 1,
                    "query_end": 68,
                    "query_start": 36,
                }
            ],
            "title": "gi|887699|gb|DQ37780 Cowpox virus 15",
        }
    ],
}

RECORD4 = {
    "query": "id4",
    "alignments": [
        {
            "length": 30000,
            "hsps": [
                {
                    "btop": "",
                    "bits": 10,
                    "sbjct_end": 1400,
                    "expect": 1e-3,
                    "sbjct": "AAAAAAAAAAA",
                    "sbjct_start": 1390,
                    "query": "TACCCTGCGGCCCGCTACGGCTGGTCTCCATTT",
                    "frame": 1,
                    "query_end": 68,
                    "query_start": 36,
                },
                {
                    "btop": "",
                    "bits": 5,
                    "sbjct_end": 1400,
                    "expect": 1e-2,
                    "sbjct": "AAAAAAAAAAA",
                    "sbjct_start": 1390,
                    "query": "TACCCTGCGGCCCGCTACGGCTGGTCTCCATTT",
                    "frame": 1,
                    "query_end": 68,
                    "query_start": 36,
                },
                {
                    "btop": "",
                    "bits": 3,
                    "sbjct_end": 1400,
                    "expect": 0.0,
                    "sbjct": "AAAAAAAAAAA",
                    "sbjct_start": 1390,
                    "query": "TACCCTGCGGCCCGCTACGGCTGGTCTCCATTT",
                    "frame": 1,
                    "query_end": 68,
                    "query_start": 36,
                },
            ],
            "title": "gi|887699|gb|DQ37780 Cowpox virus 15",
        }
    ],
}
