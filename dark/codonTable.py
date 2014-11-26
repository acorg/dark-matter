# A table with which codons translate to which amino acids.
# Based on information from: http://wang.salk.edu/research.php

CODONS = {
    'Phe': ['TTT', 'TTC'],
    'Leu': ['TTA', 'TTG', 'CTC', 'CTA'],
    'Ile': ['ATC', 'ATA'],
    'Val': ['GTC', 'GTA'],
    'Met': ['ATG'],
    'Ser': ['TCC', 'TCA', 'AGT', 'AGC'],
    'Pro': ['CCC', 'CCA'],
    'Thr': ['ACC', 'ACA'],
    'Ala': ['GCC', 'GCA'],
    'Tyr': ['TAT', 'TAC'],
    'His': ['CAT', 'CAC'],
    'Gln': ['CAA', 'CAG'],
    'Asn': ['AAT', 'AAC'],
    'Lys': ['AAA', 'AAG'],
    'Asp': ['GAT', 'GAC'],
    'Glu': ['GAA', 'GAG'],
    'Cys': ['TGT', 'TGC'],
    'Trp': ['TGG'],
    'Arg': ['CGC', 'CGA', 'AGA', 'AGG'],
    'Gly': ['GGC', 'GGA']
}
