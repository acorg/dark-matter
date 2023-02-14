# Sample BLAST parameters.
PARAMS = {
    "application": "BLASTN",
    "blast_cutoff": [None, None],
    "database": "manx-shearwater",
    "database_length": 17465129,
    "database_letters": None,
    "database_name": [],
    "database_sequences": 70016,
    "date": "",
    "dropoff_1st_pass": [None, None],
    "effective_database_length": None,
    "effective_hsp_length": 22,
    "effective_query_length": None,
    "effective_search_space": 382194648.0,
    "effective_search_space_used": None,
    "frameshift": [None, None],
    "gap_penalties": [5, 2],
    "gap_trigger": [None, None],
    "gap_x_dropoff": [None, None],
    "gap_x_dropoff_final": [None, None],
    "gapped": 0,
    "hsps_gapped": None,
    "hsps_no_gap": None,
    "hsps_prelim_gapped": None,
    "hsps_prelim_gapped_attemped": None,
    "ka_params": [0.625, 0.41, 0.78],
    "ka_params_gap": [None, None, None],
    "matrix": "",
    "num_good_extends": None,
    "num_hits": None,
    "num_letters_in_database": 17465129,
    "num_seqs_better_e": None,
    "num_sequences": None,
    "num_sequences_in_database": 70016,
    "posted_date": [],
    "query": "GZG3DGY01ASHXW",
    "query_id": "Query_1",
    "query_length": 46,
    "query_letters": 46,
    "reference": "Stephen F. Altschul, Thomas L. Madden, ...",
    "sc_match": 2,
    "sc_mismatch": -3,
    "threshold": None,
    "version": "2.2.28+",
    "window_size": None,
}

RECORD0 = {
    "query": "id0",
    "alignments": [
        {
            "length": 37000,
            "hsps": [
                {
                    "bits": 20,
                    "sbjct_end": 15400,
                    "expect": 1e-11,
                    "sbjct": "TACCC--CGGCCCGCG-CGGCCGGCTCTCCA",
                    "sbjct_start": 15362,
                    "query": "TACCCTGCGGCCCGCTACGGCTGG-TCTCCA",
                    "frame": [1, 1],
                    "query_end": 68,
                    "query_start": 28,
                }
            ],
            "title": "gi|887699|gb|DQ37780 Squirrelpox virus 1296/99",
        },
        {
            "length": 38000,
            "hsps": [
                {
                    "bits": 25,
                    "sbjct_end": 12400,
                    "expect": 1e-10,
                    "sbjct": "TACCC--CGGCCCGCG-CGGCCGGCTCTCCA",
                    "sbjct_start": 12362,
                    "query": "TACCCTGCGGCCCGCTACGGCTGG-TCTCCA",
                    "frame": [1, 1],
                    "query_end": 68,
                    "query_start": 28,
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
                    "bits": 20,
                    "sbjct_end": 11400,
                    "expect": 1e-8,
                    "sbjct": "TACCC--CGGCCCGCG-CGGCCGGCTCTCCA",
                    "sbjct_start": 11362,
                    "query": "TACCCTGCGGCCCGCTACGGCTGG-TCTCCA",
                    "frame": [1, 1],
                    "query_end": 68,
                    "query_start": 28,
                }
            ],
            "title": "gi|887699|gb|DQ37780 Monkeypox virus 456",
        },
        {
            "length": 35000,
            "hsps": [
                {
                    "bits": 20,
                    "sbjct_end": 10400,
                    "expect": 1e-7,
                    "sbjct": "TACCC--CGGCCCGCG-CGGCCGGCTCTCCA",
                    "sbjct_start": 10362,
                    "query": "TACCCTGCGGCCCGCTACGGCTGG-TCTCCA",
                    "frame": [1, 1],
                    "query_end": 68,
                    "query_start": 28,
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
                    "bits": 20,
                    "sbjct_end": 1400,
                    "expect": 1e-6,
                    "sbjct": "TACCC--CGGCCCGCG-CGGCCGGCTCTCCA",
                    "sbjct_start": 1362,
                    "query": "TACCCTGCGGCCCGCTACGGCTGG-TCTCCA",
                    "frame": [1, 1],
                    "query_end": 68,
                    "query_start": 28,
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
                    "bits": 20,
                    "sbjct_end": 1400,
                    "expect": 1e-5,
                    "sbjct": "TACCC--CGGCCCGCG-CGGCCGGCTCTCCA",
                    "sbjct_start": 1362,
                    "query": "TACCCTGCGGCCCGCTACGGCTGG-TCTCCA",
                    "frame": [1, 1],
                    "query_end": 68,
                    "query_start": 28,
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
                    "bits": 10,
                    "sbjct_end": 1400,
                    "expect": 1e-3,
                    "sbjct": "TACCC--CGGCCCGCG-CGGCCGGCTCTCCA",
                    "sbjct_start": 1362,
                    "query": "TACCCTGCGGCCCGCTACGGCTGG-TCTCCA",
                    "frame": [1, 1],
                    "query_end": 68,
                    "query_start": 28,
                },
                {
                    "bits": 5,
                    "sbjct_end": 1400,
                    "expect": 1e-2,
                    "sbjct": "TACCC--CGGCCCGCG-CGGCCGGCTCTCCA",
                    "sbjct_start": 1362,
                    "query": "TACCCTGCGGCCCGCTACGGCTGG-TCTCCA",
                    "frame": [1, 1],
                    "query_end": 68,
                    "query_start": 28,
                },
                {
                    "bits": 3,
                    "sbjct_end": 1400,
                    "expect": 0.0,
                    "sbjct": "TACCC--CGGCCCGCG-CGGCCGGCTCTCCA",
                    "sbjct_start": 1362,
                    "query": "TACCCTGCGGCCCGCTACGGCTGG-TCTCCA",
                    "frame": [1, 1],
                    "query_end": 68,
                    "query_start": 28,
                },
            ],
            "title": "gi|887699|gb|DQ37780 Cowpox virus 15",
        }
    ],
}
