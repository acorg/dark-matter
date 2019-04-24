from dark.utils import countPrint
try:
    from itertools import zip_longest
except ImportError:
    # zip_longest does not exist in Python 2.7 itertools. We should be able
    # to get it via from six.moves import zip_longest according to
    # https://pythonhosted.org/six/index.html?highlight=zip_longest but
    # that doesn't work for me.
    from itertools import izip_longest as zip_longest

# From https://en.wikipedia.org/wiki/Amino_acid
#
# Alanine          Ala     A
# Arginine         Arg     R
# Asparagine       Asn     N
# Aspartic acid    Asp     D
# Cysteine         Cys     C
# Glutamic acid    Glu     E
# Glutamine        Gln     Q
# Glycine          Gly     G
# Histidine        His     H
# Isoleucine       Ile     I
# Leucine          Leu     L
# Lysine           Lys     K
# Methionine       Met     M
# Phenylalanine    Phe     F
# Proline          Pro     P
# Serine           Ser     S
# Threonine        Thr     T
# Tryptophan       Trp     W
# Tyrosine         Tyr     Y
# Valine           Val     V

NAMES = {
    'A': 'Alanine',
    'R': 'Arginine',
    'N': 'Asparagine',
    'D': 'Aspartic acid',
    'C': 'Cysteine',
    'E': 'Glutamic acid',
    'Q': 'Glutamine',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'L': 'Leucine',
    'K': 'Lysine',
    'M': 'Methionine',
    'F': 'Phenylalanine',
    'P': 'Proline',
    'S': 'Serine',
    'T': 'Threonine',
    'V': 'Valine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine',
}

AA_LETTERS = sorted(NAMES.keys())

NAMES_TO_ABBREV1 = dict((name, abbrev1) for abbrev1, name in NAMES.items())

ABBREV3 = {
    'A': 'Ala',
    'R': 'Arg',
    'N': 'Asn',
    'D': 'Asp',
    'C': 'Cys',
    'E': 'Glu',
    'Q': 'Gln',
    'G': 'Gly',
    'H': 'His',
    'I': 'Ile',
    'L': 'Leu',
    'K': 'Lys',
    'M': 'Met',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'T': 'Thr',
    'V': 'Val',
    'W': 'Trp',
    'Y': 'Tyr',
}

ABBREV3_TO_ABBREV1 = dict((abbrev3, abbrev1)
                          for abbrev1, abbrev3 in ABBREV3.items())

HYDROPHOBIC = 0x0001
HYDROPHILIC = 0x0002
AROMATIC = 0x0004
SULPHUR = 0x0008
ALIPHATIC = 0x0010
HYDROXYLIC = 0x0020
TINY = 0x0040
SMALL = 0x0080
ACIDIC = 0x0100
BASIC_POSITIVE = 0x0200
NEGATIVE = 0x0400
POLAR = 0x0800
NONE = 0x1000


ALL_PROPERTIES = (
    ACIDIC, ALIPHATIC, AROMATIC, BASIC_POSITIVE, HYDROPHILIC,
    HYDROPHOBIC, HYDROXYLIC, NEGATIVE, NONE, POLAR, SMALL, SULPHUR, TINY)

PROPERTY_NAMES = {
    ACIDIC: 'Acidic',
    ALIPHATIC: 'Aliphatic',
    AROMATIC: 'Aromatic',
    BASIC_POSITIVE: 'Basic positive',
    HYDROPHILIC: 'Hydrophilic',
    HYDROPHOBIC: 'Hydrophobic',
    HYDROXYLIC: 'Hydroxylic',
    NEGATIVE: 'Negative',
    NONE: '<NONE>',
    POLAR: 'Polar',
    SMALL: 'Small',
    SULPHUR: 'Sulphur',
    TINY: 'Tiny',
}

PROPERTIES = {
    'A': HYDROPHOBIC | SMALL | TINY,
    'C': HYDROPHOBIC | SMALL | TINY | SULPHUR,
    'D': HYDROPHILIC | SMALL | POLAR | NEGATIVE,
    'E': HYDROPHILIC | NEGATIVE | ACIDIC,
    'F': HYDROPHOBIC | AROMATIC,
    'G': HYDROPHILIC | SMALL | TINY,
    'H': HYDROPHOBIC | AROMATIC | POLAR | BASIC_POSITIVE,
    'I': ALIPHATIC | HYDROPHOBIC,
    'K': HYDROPHOBIC | BASIC_POSITIVE | POLAR,
    'L': ALIPHATIC | HYDROPHOBIC,
    'M': HYDROPHOBIC | SULPHUR,
    'N': HYDROPHILIC | SMALL | POLAR | ACIDIC,
    'P': HYDROPHILIC | SMALL,
    'Q': HYDROPHILIC | POLAR | ACIDIC,
    'R': HYDROPHILIC | POLAR | BASIC_POSITIVE,
    'S': HYDROPHILIC | SMALL | POLAR | HYDROXYLIC,
    'T': HYDROPHOBIC | SMALL | HYDROXYLIC,
    'V': ALIPHATIC | HYDROPHOBIC | SMALL,
    'W': HYDROPHOBIC | AROMATIC | POLAR,
    'Y': HYDROPHOBIC | AROMATIC | POLAR,
}


# A table with which codons translate to which amino acids.
# Based on https://en.wikipedia.org/wiki/DNA_codon_table
#
# Note that the trailing commas are necessary in the AAs that only have one
# codon. If you omit them, the parentheses will not create a tuple.

CODONS = {
    'A': ('GCA', 'GCC', 'GCG', 'GCT',),
    'C': ('TGC', 'TGT',),
    'D': ('GAC', 'GAT',),
    'E': ('GAA', 'GAG',),
    'F': ('TTC', 'TTT',),
    'G': ('GGA', 'GGC', 'GGG', 'GGT',),
    'H': ('CAC', 'CAT',),
    'I': ('ATA', 'ATC', 'ATT',),
    'K': ('AAA', 'AAG',),
    'L': ('CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG',),
    'M': ('ATG',),
    'N': ('AAC', 'AAT',),
    'P': ('CCA', 'CCC', 'CCG', 'CCT',),
    'Q': ('CAA', 'CAG',),
    'R': ('AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT',),
    'S': ('AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT',),
    'T': ('ACA', 'ACC', 'ACG', 'ACT',),
    'V': ('GTA', 'GTC', 'GTG', 'GTT',),
    'W': ('TGG',),
    'Y': ('TAC', 'TAT',),
}

START_CODON = 'ATG'
STOP_CODONS = ('TAA', 'TAG', 'TGA',)


"""
The dictionary below contains for each amino acid the value for
each property scaled from -1 to 1.
For documentation, check https://notebooks.antigenic-cartography.org/barbara/
pages/features/aa-properties.html
"""

PROPERTY_DETAILS = {
    'A': {
        'aliphaticity': 0.305785123967,
        'aromaticity': -0.550128534704,
        'composition': -1.0,
        'hydrogenation': 0.8973042362,
        'hydropathy': 0.4,
        'hydroxythiolation': -0.265160523187,
        'iep': -0.191489361702,
        'polar requirement': -0.463414634146,
        'polarity': -0.20987654321,
        'volume': -0.664670658683,
    },
    'C': {
        'aliphaticity': -0.00826446280992,
        'aromaticity': -0.740359897172,
        'composition': 1.0,
        'hydrogenation': 0.240051347882,
        'hydropathy': 0.555555555556,
        'hydroxythiolation': 0.785969084423,
        'iep': -0.424280350438,
        'polar requirement': -1.0,
        'polarity': -0.851851851852,
        'volume': -0.377245508982,
    },
    'D': {
        'aliphaticity': -0.818181818182,
        'aromaticity': -1.0,
        'composition': 0.00363636363636,
        'hydrogenation': -0.90243902439,
        'hydropathy': -0.777777777778,
        'hydroxythiolation': -0.348394768133,
        'iep': -1.0,
        'polar requirement': 1.0,
        'polarity': 1.0,
        'volume': -0.389221556886,
    },
    'E': {
        'aliphaticity': -0.553719008264,
        'aromaticity': -0.899742930591,
        'composition': -0.330909090909,
        'hydrogenation': -1.0,
        'hydropathy': -0.777777777778,
        'hydroxythiolation': -0.555291319857,
        'iep': -0.887359198999,
        'polar requirement': 0.878048780488,
        'polarity': 0.827160493827,
        'volume': -0.0419161676647,
    },
    'F': {
        'aliphaticity': 0.223140495868,
        'aromaticity': 0.858611825193,
        'composition': -1.0,
        'hydrogenation': 0.0218228498074,
        'hydropathy': 0.622222222222,
        'hydroxythiolation': 0.0582639714625,
        'iep': -0.321652065081,
        'polar requirement': -0.951219512195,
        'polarity': -0.925925925926,
        'volume': 0.544910179641,
    },
    'G': {
        'aliphaticity': -1.0,
        'aromaticity': -0.45501285347,
        'composition': -0.461818181818,
        'hydrogenation': 1.0,
        'hydropathy': -0.0888888888889,
        'hydroxythiolation': -0.158145065398,
        'iep': -0.198998748436,
        'polar requirement': -0.243902439024,
        'polarity': 0.0123456790123,
        'volume': -1.0,
    },
    'H': {
        'aliphaticity': -0.256198347107,
        'aromaticity': 0.555269922879,
        'composition': -0.578181818182,
        'hydrogenation': -0.150192554557,
        'hydropathy': -0.711111111111,
        'hydroxythiolation': 0.0154577883472,
        'iep': 0.206508135169,
        'polar requirement': -0.121951219512,
        'polarity': 0.358024691358,
        'volume': 0.11377245509,
    },
    'I': {
        'aliphaticity': 0.867768595041,
        'aromaticity': -0.264781491003,
        'composition': -1.0,
        'hydrogenation': 0.432605905006,
        'hydropathy': 1.0,
        'hydroxythiolation': -0.85255648038,
        'iep': -0.18648310388,
        'polar requirement': -0.975609756098,
        'polarity': -0.925925925926,
        'volume': 0.293413173653,
    },
    'K': {
        'aliphaticity': 0.123966942149,
        'aromaticity': -0.141388174807,
        'composition': -0.76,
        'hydrogenation': -0.142490372272,
        'hydropathy': -0.866666666667,
        'hydroxythiolation': -1.0,
        'iep': 0.744680851064,
        'polar requirement': 0.292682926829,
        'polarity': 0.58024691358,
        'volume': 0.389221556886,
    },
    'L': {
        'aliphaticity': 1.0,
        'aromaticity': -0.287917737789,
        'composition': -1.0,
        'hydrogenation': 0.381258023107,
        'hydropathy': 0.844444444444,
        'hydroxythiolation': -0.745541022592,
        'iep': -0.196495619524,
        'polar requirement': -0.975609756098,
        'polarity': -1.0,
        'volume': 0.293413173653,
    },
    'M': {
        'aliphaticity': 0.537190082645,
        'aromaticity': -0.372750642674,
        'composition': -1.0,
        'hydrogenation': -0.186136071887,
        'hydropathy': 0.422222222222,
        'hydroxythiolation': 0.0653983353151,
        'iep': -0.256570713392,
        'polar requirement': -0.878048780488,
        'polarity': -0.802469135802,
        'volume': 0.221556886228,
    },
    'N': {
        'aliphaticity': 0.471074380165,
        'aromaticity': -0.616966580977,
        'composition': -0.0327272727273,
        'hydrogenation': -0.548138639281,
        'hydropathy': -0.777777777778,
        'hydroxythiolation': 0.277051129608,
        'iep': -0.339173967459,
        'polar requirement': 0.268292682927,
        'polarity': 0.654320987654,
        'volume': -0.365269461078,
    },
    'P': {
        'aliphaticity': -0.917355371901,
        'aromaticity': -0.308483290488,
        'composition': -0.716363636364,
        'hydrogenation': 1.0,
        'hydropathy': -0.355555555556,
        'hydroxythiolation': -0.203329369798,
        'iep': -0.116395494368,
        'polar requirement': -0.560975609756,
        'polarity': -0.234567901235,
        'volume': -0.646706586826,
    },
    'Q': {
        'aliphaticity': 0.652892561983,
        'aromaticity': -0.439588688946,
        'composition': -0.352727272727,
        'hydrogenation': -0.602053915276,
        'hydropathy': -0.777777777778,
        'hydroxythiolation': -0.177170035672,
        'iep': -0.279098873592,
        'polar requirement': -0.0731707317073,
        'polarity': 0.382716049383,
        'volume': -0.0179640718563,
    },
    'R': {
        'aliphaticity': -0.157024793388,
        'aromaticity': -0.0642673521851,
        'composition': -0.527272727273,
        'hydrogenation': -0.401797175866,
        'hydropathy': -1.0,
        'hydroxythiolation': -0.51486325802,
        'iep': 1.0,
        'polar requirement': 0.0487804878049,
        'polarity': 0.382716049383,
        'volume': 0.449101796407,
    },
    'S': {
        'aliphaticity': 0.256198347107,
        'aromaticity': -0.660668380463,
        'composition': 0.0327272727273,
        'hydrogenation': 0.106546854942,
        'hydropathy': -0.177777777778,
        'hydroxythiolation': 1.0,
        'iep': -0.271589486859,
        'polar requirement': -0.341463414634,
        'polarity': 0.0617283950617,
        'volume': -0.652694610778,
    },
    'T': {
        'aliphaticity': -0.123966942149,
        'aromaticity': -0.80205655527,
        'composition': -0.483636363636,
        'hydrogenation': 0.399229781772,
        'hydropathy': -0.155555555556,
        'hydroxythiolation': 0.709869203329,
        'iep': -0.151439299124,
        'polar requirement': -0.560975609756,
        'polarity': -0.0864197530864,
        'volume': -0.305389221557,
    },
    'V': {
        'aliphaticity': 0.570247933884,
        'aromaticity': -0.665809768638,
        'composition': -1.0,
        'hydrogenation': 0.679075738126,
        'hydropathy': 0.933333333333,
        'hydroxythiolation': -0.621878715815,
        'iep': -0.201501877347,
        'polar requirement': -0.80487804878,
        'polarity': -0.753086419753,
        'volume': -0.0299401197605,
    },
    'W': {
        'aliphaticity': -0.619834710744,
        'aromaticity': 1.0,
        'composition': -0.905454545455,
        'hydrogenation': 0.0218228498074,
        'hydropathy': -0.2,
        'hydroxythiolation': 0.00118906064209,
        'iep': -0.219023779725,
        'polar requirement': -0.90243902439,
        'polarity': -0.876543209877,
        'volume': 1.0,
    },
    'Y': {
        'aliphaticity': -0.454545454545,
        'aromaticity': 0.712082262211,
        'composition': -0.854545454545,
        'hydrogenation': -0.304236200257,
        'hydropathy': 0.288888888889,
        'hydroxythiolation': 0.405469678954,
        'iep': -0.276595744681,
        'polar requirement': -0.853658536585,
        'polarity': -0.679012345679,
        'volume': 0.592814371257,
    },
}


"""
The dictionary below contains for each amino acid the value for
each property.
"""

PROPERTY_DETAILS_RAW = {
    'A': {
        'aliphaticity': 0.239,
        'aromaticity': -0.11,
        'composition': 0.0,
        'hydrogenation': 0.33,
        'hydropathy': 1.8,
        'hydroxythiolation': -0.062,
        'iep': 6.0,
        'polar requirement': 7.0,
        'polarity': 8.1,
        'volume': 31.0,
    },
    'C': {
        'aliphaticity': 0.22,
        'aromaticity': -0.184,
        'composition': 2.75,
        'hydrogenation': 0.074,
        'hydropathy': 2.5,
        'hydroxythiolation': 0.38,
        'iep': 5.07,
        'polar requirement': 4.8,
        'polarity': 5.5,
        'volume': 55.0,
    },
    'D': {
        'aliphaticity': 0.171,
        'aromaticity': -0.285,
        'composition': 1.38,
        'hydrogenation': -0.371,
        'hydropathy': -3.5,
        'hydroxythiolation': -0.079,
        'iep': 2.77,
        'polar requirement': 13.0,
        'polarity': 13.0,
        'volume': 54.0,
    },
    'E': {
        'aliphaticity': 0.187,
        'aromaticity': -0.246,
        'composition': 0.92,
        'hydrogenation': -0.409,
        'hydropathy': -3.5,
        'hydroxythiolation': -0.184,
        'iep': 3.22,
        'polar requirement': 12.5,
        'polarity': 12.3,
        'volume': 83.0,
    },
    'F': {
        'aliphaticity': 0.234,
        'aromaticity': 0.438,
        'composition': 0.0,
        'hydrogenation': -0.011,
        'hydropathy': 2.8,
        'hydroxythiolation': 0.074,
        'iep': 5.48,
        'polar requirement': 5.0,
        'polarity': 5.4,
        'volume': 132.0,
    },
    'G': {
        'aliphaticity': 0.16,
        'aromaticity': -0.073,
        'composition': 0.74,
        'hydrogenation': 0.37,
        'hydropathy': -0.4,
        'hydroxythiolation': -0.017,
        'iep': 5.97,
        'polar requirement': 7.9,
        'polarity': 9.0,
        'volume': 3.0,
    },
    'H': {
        'aliphaticity': 0.205,
        'aromaticity': 0.32,
        'composition': 0.58,
        'hydrogenation': -0.078,
        'hydropathy': -3.2,
        'hydroxythiolation': 0.056,
        'iep': 7.59,
        'polar requirement': 8.4,
        'polarity': 10.4,
        'volume': 96.0,
    },
    'I': {
        'aliphaticity': 0.273,
        'aromaticity': 0.001,
        'composition': 0.0,
        'hydrogenation': 0.149,
        'hydropathy': 4.5,
        'hydroxythiolation': -0.309,
        'iep': 6.02,
        'polar requirement': 4.9,
        'polarity': 5.2,
        'volume': 111.0,
    },
    'K': {
        'aliphaticity': 0.228,
        'aromaticity': 0.049,
        'composition': 0.33,
        'hydrogenation': -0.075,
        'hydropathy': -3.9,
        'hydroxythiolation': -0.371,
        'iep': 9.74,
        'polar requirement': 10.1,
        'polarity': 11.3,
        'volume': 119.0,
    },
    'L': {
        'aliphaticity': 0.281,
        'aromaticity': -0.008,
        'composition': 0.0,
        'hydrogenation': 0.129,
        'hydropathy': 3.8,
        'hydroxythiolation': -0.264,
        'iep': 5.98,
        'polar requirement': 4.9,
        'polarity': 4.9,
        'volume': 111.0,
    },
    'M': {
        'aliphaticity': 0.253,
        'aromaticity': -0.041,
        'composition': 0.0,
        'hydrogenation': -0.092,
        'hydropathy': 1.9,
        'hydroxythiolation': 0.077,
        'iep': 5.74,
        'polar requirement': 5.3,
        'polarity': 5.7,
        'volume': 105.0,
    },
    'N': {
        'aliphaticity': 0.249,
        'aromaticity': -0.136,
        'composition': 1.33,
        'hydrogenation': -0.233,
        'hydropathy': -3.5,
        'hydroxythiolation': 0.166,
        'iep': 5.41,
        'polar requirement': 10.0,
        'polarity': 11.6,
        'volume': 56.0,
    },
    'P': {
        'aliphaticity': 0.165,
        'aromaticity': -0.016,
        'composition': 0.39,
        'hydrogenation': 0.37,
        'hydropathy': -1.6,
        'hydroxythiolation': -0.036,
        'iep': 6.3,
        'polar requirement': 6.6,
        'polarity': 8.0,
        'volume': 32.5,
    },
    'Q': {
        'aliphaticity': 0.26,
        'aromaticity': -0.067,
        'composition': 0.89,
        'hydrogenation': -0.254,
        'hydropathy': -3.5,
        'hydroxythiolation': -0.025,
        'iep': 5.65,
        'polar requirement': 8.6,
        'polarity': 10.5,
        'volume': 85.0,
    },
    'R': {
        'aliphaticity': 0.211,
        'aromaticity': 0.079,
        'composition': 0.65,
        'hydrogenation': -0.176,
        'hydropathy': -4.5,
        'hydroxythiolation': -0.167,
        'iep': 10.76,
        'polar requirement': 9.1,
        'polarity': 10.5,
        'volume': 124.0,
    },
    'S': {
        'aliphaticity': 0.236,
        'aromaticity': -0.153,
        'composition': 1.42,
        'hydrogenation': 0.022,
        'hydropathy': -0.8,
        'hydroxythiolation': 0.47,
        'iep': 5.68,
        'polar requirement': 7.5,
        'polarity': 9.2,
        'volume': 32.0,
    },
    'T': {
        'aliphaticity': 0.213,
        'aromaticity': -0.208,
        'composition': 0.71,
        'hydrogenation': 0.136,
        'hydropathy': -1.3,
        'hydroxythiolation': 0.348,
        'iep': 6.16,
        'polar requirement': 6.6,
        'polarity': 8.6,
        'volume': 61.0,
    },
    'V': {
        'aliphaticity': 0.255,
        'aromaticity': -0.155,
        'composition': 0.0,
        'hydrogenation': 0.245,
        'hydropathy': 4.2,
        'hydroxythiolation': -0.212,
        'iep': 5.96,
        'polar requirement': 5.6,
        'polarity': 5.9,
        'volume': 84.0,
    },
    'W': {
        'aliphaticity': 0.183,
        'aromaticity': 0.493,
        'composition': 0.13,
        'hydrogenation': -0.011,
        'hydropathy': -0.9,
        'hydroxythiolation': 0.05,
        'iep': 5.89,
        'polar requirement': 5.2,
        'polarity': 5.4,
        'volume': 170.0,
    },
    'Y': {
        'aliphaticity': 0.193,
        'aromaticity': 0.183,
        'composition': 0.2,
        'hydrogenation': -0.138,
        'hydropathy': -1.3,
        'hydroxythiolation': 0.22,
        'iep': 5.66,
        'polar requirement': 5.4,
        'polarity': 6.2,
        'volume': 136.0,
    },
}


"""
Clusters based on raw amino acid property values. See
https://notebooks.antigenic-cartography.org/barbara/pages/features/
aa-properties.html and https://notebooks.antigenic-cartography.org/barbara/
pages/features/new-tps.html
"""

PROPERTY_CLUSTERS = {
    'A': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 3,
        'hydroxythiolation': 2,
        'iep': 2,
        'polar requirement': 2,
        'polarity': 2,
        'volume': 2,
    },
    'C': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 3,
        'hydrogenation': 1,
        'hydropathy': 3,
        'hydroxythiolation': 5,
        'iep': 2,
        'polar requirement': 1,
        'polarity': 1,
        'volume': 3,
    },
    'D': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 2,
        'hydrogenation': 1,
        'hydropathy': 1,
        'hydroxythiolation': 2,
        'iep': 1,
        'polar requirement': 4,
        'polarity': 4,
        'volume': 3,
    },
    'E': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 1,
        'hydroxythiolation': 1,
        'iep': 1,
        'polar requirement': 4,
        'polarity': 4,
        'volume': 4,
    },
    'F': {
        'aliphaticity': 1,
        'aromaticity': 2,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 3,
        'hydroxythiolation': 3,
        'iep': 2,
        'polar requirement': 1,
        'polarity': 1,
        'volume': 4,
    },
    'G': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 2,
        'hydroxythiolation': 2,
        'iep': 2,
        'polar requirement': 2,
        'polarity': 2,
        'volume': 1,
    },
    'H': {
        'aliphaticity': 1,
        'aromaticity': 2,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 1,
        'hydroxythiolation': 3,
        'iep': 3,
        'polar requirement': 2,
        'polarity': 3,
        'volume': 4,
    },
    'I': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 4,
        'hydroxythiolation': 1,
        'iep': 2,
        'polar requirement': 1,
        'polarity': 1,
        'volume': 4,
    },
    'K': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 1,
        'hydroxythiolation': 1,
        'iep': 3,
        'polar requirement': 3,
        'polarity': 4,
        'volume': 4,
    },
    'L': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 4,
        'hydroxythiolation': 1,
        'iep': 2,
        'polar requirement': 1,
        'polarity': 1,
        'volume': 4,
    },
    'M': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 3,
        'hydroxythiolation': 3,
        'iep': 2,
        'polar requirement': 1,
        'polarity': 1,
        'volume': 4,
    },
    'N': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 2,
        'hydrogenation': 1,
        'hydropathy': 1,
        'hydroxythiolation': 4,
        'iep': 2,
        'polar requirement': 3,
        'polarity': 4,
        'volume': 3,
    },
    'P': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 2,
        'hydroxythiolation': 2,
        'iep': 2,
        'polar requirement': 2,
        'polarity': 2,
        'volume': 2,
    },
    'Q': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 1,
        'hydroxythiolation': 2,
        'iep': 2,
        'polar requirement': 2,
        'polarity': 3,
        'volume': 4,
    },
    'R': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 1,
        'hydroxythiolation': 1,
        'iep': 3,
        'polar requirement': 2,
        'polarity': 3,
        'volume': 4,
    },
    'S': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 2,
        'hydrogenation': 1,
        'hydropathy': 2,
        'hydroxythiolation': 5,
        'iep': 2,
        'polar requirement': 2,
        'polarity': 2,
        'volume': 2,
    },
    'T': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 2,
        'hydroxythiolation': 5,
        'iep': 2,
        'polar requirement': 2,
        'polarity': 2,
        'volume': 3,
    },
    'V': {
        'aliphaticity': 1,
        'aromaticity': 1,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 4,
        'hydroxythiolation': 1,
        'iep': 2,
        'polar requirement': 1,
        'polarity': 1,
        'volume': 4,
    },
    'W': {
        'aliphaticity': 1,
        'aromaticity': 2,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 2,
        'hydroxythiolation': 3,
        'iep': 2,
        'polar requirement': 1,
        'polarity': 1,
        'volume': 5,
    },
    'Y': {
        'aliphaticity': 1,
        'aromaticity': 2,
        'composition': 1,
        'hydrogenation': 1,
        'hydropathy': 2,
        'hydroxythiolation': 4,
        'iep': 2,
        'polar requirement': 1,
        'polarity': 1,
        'volume': 4,
    },
}


class AminoAcid(object):
    """
    Hold information about an amino acid.

    @param name: The full C{str} name of the amino acid.
    @param abbrev3: The 3-letter C{str} abbreviation of the amino acid,
        e.g., 'Arg'.
    @param abbrev1: The 1-letter C{str} abbreviation of the amino acid,
        e.g., 'A'.
    @param codons: A C{list} of 3-letter codons for the amino acid.
    @param properties: An C{int} logical-AND of the various properties
        (see PROPERTIES, above) of this amino acid.
    @param propertyDetails: A C{dict} containing property names and values
        for this amino acid. E.g.:
        {
            'aliphaticity': -0.157024793388,
            'aromaticity': -0.0642673521851,
            'composition': -0.527272727273,
            'hydrogenation': -0.401797175866,
            'hydropathy': -1.0,
            'hydroxythiolation': -0.51486325802,
            'iep': 1.0,
            'polar requirement': 0.0487804878049,
            'polarity': 0.382716049383,
            'volume': 0.449101796407,
        }
    @param propertyClusters: A C{dict} containing the property names and
        clusters for this amino acid. E.g.:
        {
            'aliphaticity': 1,
            'aromaticity': 2,
            'composition': 1,
            'hydrogenation': 1,
            'hydropathy': 2,
            'hydroxythiolation': 4,
            'iep': 2,
            'polar requirement': 1,
            'polarity': 1,
            'volume': 4,
        }
    """
    def __init__(self, name, abbrev3, abbrev1, codons, properties,
                 propertyDetails, propertyClusters):
        self.name = name
        self.abbrev3 = abbrev3
        self.abbrev1 = abbrev1
        self.codons = codons
        self.properties = properties
        self.propertyDetails = propertyDetails
        self.propertyClusters = propertyClusters


def find(s):
    """
    Find an amino acid whose name or abbreviation is s.

    @param s: A C{str} amino acid specifier. This may be a full name,
        a 3-letter abbreviation or a 1-letter abbreviation. Case is ignored.
    return: An L{AminoAcid} instance or C{None} if no matching amino acid can
        be located.
    """

    abbrev1 = None
    origS = s

    if ' ' in s:
        # Convert first word to title case, others to lower.
        first, rest = s.split(' ', 1)
        s = first.title() + ' ' + rest.lower()
    else:
        s = s.title()

    if s in NAMES:
        abbrev1 = s
    elif s in ABBREV3_TO_ABBREV1:
        abbrev1 = ABBREV3_TO_ABBREV1[s]
    elif s in NAMES_TO_ABBREV1:
        abbrev1 = NAMES_TO_ABBREV1[s]
    else:
        # Look for a 3-letter codon.
        def findCodon(target):
            for abbrev1, codons in CODONS.items():
                for codon in codons:
                    if codon == target:
                        return abbrev1

        abbrev1 = findCodon(origS.upper())

    if abbrev1:
        return AminoAcid(
            NAMES[abbrev1], ABBREV3[abbrev1], abbrev1, CODONS[abbrev1],
            PROPERTIES[abbrev1], PROPERTY_DETAILS[abbrev1],
            PROPERTY_CLUSTERS[abbrev1])


def _propertiesOrClustersForSequence(sequence, propertyNames, propertyValues,
                                     missingAAValue):
    """
    Extract amino acid property values or cluster numbers for a sequence.

    @param sequence: An C{AARead} (or a subclass) instance.
    @param propertyNames: An iterable of C{str} property names (each of which
        must be a key of a key in the C{propertyValues} C{dict}).
    @param propertyValues: A C{dict} in the form of C{PROPERTY_DETAILS} or
        C{PROPERTY_CLUSTERS} (see above).
    @param missingAAValue: A C{float} value to use for properties when an AA
        (e.g., 'X') is not known.
    @raise ValueError: If an unknown property is given in C{propertyNames}.
    @return: A C{dict} keyed by (lowercase) property name, with values that are
        C{list}s of the corresponding property value in C{propertyValues} in
        order of sequence position.
    """
    propertyNames = sorted(map(str.lower, set(propertyNames)))

    # Make sure all mentioned property names exist for at least one AA.
    knownProperties = set()
    for names in propertyValues.values():
        knownProperties.update(names)
    unknown = set(propertyNames) - knownProperties
    if unknown:
        raise ValueError(
            'Unknown propert%s: %s.' %
            ('y' if len(unknown) == 1 else 'ies', ', '.join(unknown)))

    aas = sequence.sequence.upper()
    result = {}

    for propertyName in propertyNames:
        result[propertyName] = values = []
        append = values.append
        for aa in aas:
            try:
                properties = propertyValues[aa]
            except KeyError:
                # No such AA.
                append(missingAAValue)
            else:
                append(properties[propertyName])

    return result


def propertiesForSequence(sequence, propertyNames, missingAAValue=-1.1):
    """
    Extract amino acid property values for a sequence.

    @param sequence: An C{AARead} (or a subclass) instance.
    @param propertyNames: An iterable of C{str} property names (each of which
        must be a key of a key in the C{dark.aa.PROPERTY_DETAILS} C{dict}).
    @param missingAAValue: A C{float} value to use for properties when an AA
        (e.g., 'X') is not known.
    @raise ValueError: If an unknown property is given in C{propertyNames}.
    @return: A C{dict} keyed by (lowercase) property name, with values that are
        C{list}s of the corresponding property value according to sequence
        position.
    """
    return _propertiesOrClustersForSequence(
        sequence, propertyNames, PROPERTY_DETAILS, missingAAValue)


def clustersForSequence(sequence, propertyNames, missingAAValue=0):
    """
    Extract amino acid property cluster numbers for a sequence.

    @param sequence: An C{AARead} (or a subclass) instance.
    @param propertyNames: An iterable of C{str} property names (each of which
        must be a key of a key in the C{dark.aa.PROPERTY_CLUSTERS} C{dict}).
    @param missingAAValue: An C{int} value to use for properties when an AA
        (e.g., 'X') is not known.
    @raise ValueError: If an unknown property is given in C{propertyNames}.
    @return: A C{dict} keyed by (lowercase) property name, with values that are
        C{list}s of the corresponding property cluster number according to
        sequence position.
    """
    return _propertiesOrClustersForSequence(
        sequence, propertyNames, PROPERTY_CLUSTERS, missingAAValue)


def matchToString(aaMatch, read1, read2, indent='', offsets=None):
    """
    Format amino acid sequence match as a string.

    @param aaMatch: A C{dict} returned by C{compareAaReads}.
    @param read1: A C{Read} instance or an instance of one of its subclasses.
    @param read2: A C{Read} instance or an instance of one of its subclasses.
    @param indent: A C{str} to indent all returned lines with.
    @param offsets: If not C{None}, a C{set} of offsets of interest that were
        only considered when making C{match}.
    @return: A C{str} describing the match.
    """
    match = aaMatch['match']
    matchCount = match['matchCount']
    gapMismatchCount = match['gapMismatchCount']
    gapGapMismatchCount = match['gapGapMismatchCount']
    nonGapMismatchCount = match['nonGapMismatchCount']

    if offsets:
        len1 = len2 = len(offsets)
    else:
        len1, len2 = map(len, (read1, read2))

    result = []
    append = result.append

    append(countPrint('%sMatches' % indent, matchCount, len1, len2))
    mismatchCount = (gapMismatchCount + gapGapMismatchCount +
                     nonGapMismatchCount)
    append(countPrint('%sMismatches' % indent, mismatchCount, len1, len2))
    append(countPrint('%s  Not involving gaps (i.e., conflicts)' % (indent),
                      nonGapMismatchCount, len1, len2))
    append(countPrint('%s  Involving a gap in one sequence' % indent,
                      gapMismatchCount, len1, len2))
    append(countPrint('%s  Involving a gap in both sequences' % indent,
                      gapGapMismatchCount, len1, len2))

    for read, key in zip((read1, read2), ('read1', 'read2')):
        append('%s  Id: %s' % (indent, read.id))
        length = len(read)
        append('%s    Length: %d' % (indent, length))
        gapCount = len(aaMatch[key]['gapOffsets'])
        append(countPrint('%s    Gaps' % indent, gapCount, length))
        if gapCount:
            append(
                '%s    Gap locations (1-based): %s' %
                (indent,
                 ', '.join(map(lambda offset: str(offset + 1),
                               sorted(aaMatch[key]['gapOffsets'])))))
        extraCount = aaMatch[key]['extraCount']
        if extraCount:
            append(countPrint('%s    Extra nucleotides at end' % indent,
                              extraCount, length))

    return '\n'.join(result)


def compareAaReads(read1, read2, gapChars='-', offsets=None):
    """
    Compare two amino acid sequences.

    @param read1: A C{Read} instance or an instance of one of its subclasses.
    @param read2: A C{Read} instance or an instance of one of its subclasses.
    @param gapChars: An object supporting __contains__ with characters that
        should be considered to be gaps.
    @param offsets: If not C{None}, a C{set} of offsets of interest. Offsets
        not in the set will not be considered.
    @return: A C{dict} with information about the match and the individual
        sequences (see below).
    """
    matchCount = 0
    gapMismatchCount = nonGapMismatchCount = gapGapMismatchCount = 0
    read1ExtraCount = read2ExtraCount = 0
    read1GapOffsets = []
    read2GapOffsets = []

    for offset, (a, b) in enumerate(zip_longest(read1.sequence.upper(),
                                                read2.sequence.upper())):
        # Use 'is not None' in the following to allow an empty offsets set
        # to be passed.
        if offsets is not None and offset not in offsets:
            continue
        if a is None:
            # b has an extra character at its end (it cannot be None).
            assert b is not None
            read2ExtraCount += 1
            if b in gapChars:
                read2GapOffsets.append(offset)
        elif b is None:
            # a has an extra character at its end.
            read1ExtraCount += 1
            if a in gapChars:
                read1GapOffsets.append(offset)
        else:
            # We have a character from both sequences (they could still be
            # gap characters).
            if a in gapChars:
                read1GapOffsets.append(offset)
                if b in gapChars:
                    # Both are gaps. This can happen (though hopefully not
                    # if the sequences were pairwise aligned).
                    gapGapMismatchCount += 1
                    read2GapOffsets.append(offset)
                else:
                    # a is a gap, b is not.
                    gapMismatchCount += 1
            else:
                if b in gapChars:
                    # b is a gap, a is not.
                    gapMismatchCount += 1
                    read2GapOffsets.append(offset)
                else:
                    # Neither is a gap character.
                    if a == b:
                        matchCount += 1
                    else:
                        nonGapMismatchCount += 1

    return {
        'match': {
            'matchCount': matchCount,
            'gapMismatchCount': gapMismatchCount,
            'gapGapMismatchCount': gapGapMismatchCount,
            'nonGapMismatchCount': nonGapMismatchCount,
        },
        'read1': {
            'extraCount': read1ExtraCount,
            'gapOffsets': read1GapOffsets,
        },
        'read2': {
            'extraCount': read2ExtraCount,
            'gapOffsets': read2GapOffsets,
        },
    }
