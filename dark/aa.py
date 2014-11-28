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

PROPERTIES = {
    'I': ALIPHATIC | HYDROPHOBIC,
    'L': ALIPHATIC | HYDROPHOBIC,
    'V': ALIPHATIC | HYDROPHOBIC | SMALL,
    'M': HYDROPHOBIC | SULPHUR,
    'F': HYDROPHOBIC | AROMATIC,
    'A': HYDROPHOBIC | SMALL | TINY,
    'C': HYDROPHOBIC | SMALL | TINY | SULPHUR,
    'T': HYDROPHOBIC | SMALL | HYDROXYLIC,
    'Y': HYDROPHOBIC | AROMATIC | POLAR,
    'W': HYDROPHOBIC | AROMATIC | POLAR,
    'H': HYDROPHOBIC | AROMATIC | POLAR | BASIC_POSITIVE,
    'K': HYDROPHOBIC | BASIC_POSITIVE | POLAR,
    'P': HYDROPHILIC | SMALL,
    'G': HYDROPHILIC | SMALL | TINY,
    'S': HYDROPHILIC | SMALL | POLAR | HYDROXYLIC,
    'N': HYDROPHILIC | SMALL | POLAR | ACIDIC,
    'D': HYDROPHILIC | SMALL | POLAR | NEGATIVE,
    'Q': HYDROPHILIC | POLAR | ACIDIC,
    'E': HYDROPHILIC | NEGATIVE | ACIDIC,
    'R': HYDROPHILIC | POLAR | BASIC_POSITIVE
}


# A table with which codons translate to which amino acids.
# Based on information from: http://wang.salk.edu/research.php

CODONS = {
    'phe': ['TTT', 'TTC'],
    'leu': ['TTA', 'TTG', 'CTC', 'CTA'],
    'ile': ['ATC', 'ATA'],
    'val': ['GTC', 'GTA'],
    'met': ['ATG'],
    'ser': ['TCC', 'TCA', 'AGT', 'AGC'],
    'pro': ['CCC', 'CCA'],
    'thr': ['ACC', 'ACA'],
    'ala': ['GCC', 'GCA'],
    'tyr': ['TAT', 'TAC'],
    'his': ['CAT', 'CAC'],
    'gln': ['CAA', 'CAG'],
    'asn': ['AAT', 'AAC'],
    'lys': ['AAA', 'AAG'],
    'asp': ['GAT', 'GAC'],
    'glu': ['GAA', 'GAG'],
    'cys': ['TGT', 'TGC'],
    'trp': ['TGG'],
    'arg': ['CGC', 'CGA', 'AGA', 'AGG'],
    'gly': ['GGC', 'GGA'],
}


"""
The dictionary below contains for each amino acid the value for
each property scaled from -1 to 1.
"""

PROPERTY_DETAILS = {
    'I': {'volume': 0.293413173653, 'iep': -0.18648310388,
          'polarity': -0.925925925926,
          'composition': -1.0, 'polar_req': -0.975609756098, 'hydropathy': 1.0,
          'aliphaticity': 0.867768595041, 'hydrogenation': 0.432605905006,
          'aromaticity': -0.264781491003, 'hydroxyethilation': -0.85255648038},
    'L': {'volume': 0.293413173653, 'iep': -0.196495619524, 'polarity': -1.0,
          'composition': -1.0, 'polar_req': -0.975609756098,
          'hydropathy': 0.844444444444, 'aliphaticity': 1.0,
          'hydrogenation': 0.381258023107, 'aromaticity': -0.287917737789,
          'hydroxyethilation': -0.745541022592},
    'V': {'volume': -0.0299401197605, 'iep': -0.201501877347,
          'polarity': -0.753086419753, 'composition': -1.0,
          'polar_req': -0.80487804878, 'hydropathy': 0.933333333333,
          'aliphaticity': 0.570247933884, 'hydrogenation': 0.679075738126,
          'aromaticity': -0.665809768638,
          'hydroxyethilation': -0.621878715815},
    'M': {'volume': 0.221556886228, 'iep': -0.256570713392, 'polarity': -0.802469135802,
          'composition': -1.0, 'polar_req': -0.878048780488, 'hydropathy': 0.422222222222,
          'aliphaticity': 0.537190082645, 'hydrogenation': -0.186136071887, 'aromaticity': -0.372750642674,
          'hydroxyethilation': 0.0653983353151},
    'F': {'volume': 0.544910179641, 'iep': -0.321652065081, 'polarity': -0.925925925926,
          'composition': -1.0, 'polar_req': -0.951219512195, 'hydropathy': 0.622222222222,
          'aliphaticity': 0.223140495868, 'hydrogenation': 0.0218228498074, 'aromaticity': 0.858611825193,
          'hydroxyethilation': 0.0582639714625},
    'A': {'volume': -0.664670658683, 'iep': -0.191489361702, 'polarity': -0.20987654321,
          'composition': -1.0, 'polar_req': -0.463414634146, 'hydropathy': 0.4,
          'aliphaticity': 0.305785123967, 'hydrogenation': 0.8973042362, 'aromaticity': -0.550128534704,
          'hydroxyethilation': -0.265160523187},
    'C': {'volume': -0.377245508982, 'iep': -0.424280350438, 'polarity': -0.851851851852,
          'composition': 1.0, 'polar_req': -1.0, 'hydropathy': 0.555555555556,
          'aliphaticity': -0.00826446280992, 'hydrogenation': 0.240051347882, 'aromaticity': -0.740359897172,
          'hydroxyethilation': 0.785969084423},
    'T': {'volume': -0.305389221557, 'iep': -0.151439299124, 'polarity': -0.0864197530864,
          'composition': -0.483636363636, 'polar_req': -0.560975609756, 'hydropathy': -0.155555555556,
          'aliphaticity': -0.123966942149, 'hydrogenation': 0.399229781772, 'aromaticity': -0.80205655527,
          'hydroxyethilation': 0.709869203329},
    'Y': {'volume': 0.592814371257, 'iep': -0.276595744681, 'polarity': -0.679012345679,
          'composition': -0.854545454545, 'polar_req': -0.853658536585, 'hydropathy': 0.288888888889,
          'aliphaticity': -0.454545454545, 'hydrogenation': -0.304236200257, 'aromaticity': 0.712082262211,
          'hydroxyethilation': 0.405469678954},
    'W': {'volume': 1.0, 'iep': -0.219023779725, 'polarity': -0.876543209877,
          'composition': -0.905454545455, 'polar_req': -0.90243902439, 'hydropathy': -0.2,
          'aliphaticity': -0.619834710744, 'hydrogenation': 0.0218228498074, 'aromaticity': 1.0,
          'hydroxyethilation': 0.00118906064209},
    'H': {'volume': 0.11377245509, 'iep': 0.206508135169, 'polarity': 0.358024691358,
          'composition': -0.578181818182, 'polar_req': -0.121951219512, 'hydropathy': -0.711111111111,
          'aliphaticity': -0.256198347107, 'hydrogenation': -0.150192554557, 'aromaticity': 0.555269922879,
          'hydroxyethilation': 0.0154577883472},
    'K': {'volume': 0.389221556886, 'iep': 0.744680851064, 'polarity': 0.58024691358,
          'composition': -0.76, 'polar_req': 0.292682926829, 'hydropathy': -0.866666666667,
          'aliphaticity': 0.123966942149, 'hydrogenation': -0.142490372272, 'aromaticity': -0.141388174807,
          'hydroxyethilation': -1.0},
    'P': {'volume': -0.646706586826, 'iep': -0.116395494368, 'polarity': -0.234567901235,
          'composition': -0.716363636364, 'polar_req': -0.560975609756, 'hydropathy': -0.355555555556,
          'aliphaticity': -0.917355371901, 'hydrogenation': 1.0, 'aromaticity': -0.308483290488,
          'hydroxyethilation': -0.203329369798},
    'G': {'volume': -1.0, 'iep': -0.198998748436, 'polarity': 0.0123456790123,
          'composition': -0.461818181818, 'polar_req': -0.243902439024, 'hydropathy': -0.0888888888889,
          'aliphaticity': -1.0, 'hydrogenation': 1.0, 'aromaticity': -0.45501285347,
          'hydroxyethilation': -0.158145065398},
    'S': {'volume': -0.652694610778, 'iep': -0.271589486859, 'polarity': 0.0617283950617,
          'composition': 0.0327272727273, 'polar_req': -0.341463414634, 'hydropathy': -0.177777777778,
          'aliphaticity': 0.256198347107, 'hydrogenation': 0.106546854942, 'aromaticity': -0.660668380463,
          'hydroxyethilation': 1.0},
    'N': {'volume': -0.365269461078, 'iep': -0.339173967459, 'polarity': 0.654320987654,
          'composition': -0.0327272727273, 'polar_req': 0.268292682927, 'hydropathy': -0.777777777778,
          'aliphaticity': 0.471074380165, 'hydrogenation': -0.548138639281, 'aromaticity': -0.616966580977,
          'hydroxyethilation': 0.277051129608},
    'D': {'volume': -0.389221556886, 'iep': -1.0, 'polarity': 1.0,
          'composition': 0.00363636363636, 'polar_req': 1.0, 'hydropathy': -0.777777777778,
          'aliphaticity': -0.818181818182, 'hydrogenation': -0.90243902439, 'aromaticity': -1.0,
          'hydroxyethilation': -0.348394768133},
    'Q': {'volume': -0.0179640718563, 'iep': -0.279098873592, 'polarity': 0.382716049383,
          'composition': -0.352727272727, 'polar_req': -0.0731707317073, 'hydropathy': -0.777777777778,
          'aliphaticity': 0.652892561983, 'hydrogenation': -0.602053915276, 'aromaticity': -0.439588688946,
          'hydroxyethilation': -0.177170035672},
    'E': {'volume': -0.0419161676647, 'iep': -0.887359198999, 'polarity': 0.827160493827,
          'composition': -0.330909090909, 'polar_req': 0.878048780488, 'hydropathy': -0.777777777778,
          'aliphaticity': -0.553719008264, 'hydrogenation': -1.0, 'aromaticity': -0.899742930591,
          'hydroxyethilation': -0.555291319857},
    'R': {'volume': 0.449101796407, 'iep': 1.0, 'polarity': 0.382716049383,
          'composition': -0.527272727273, 'polar_req': 0.0487804878049, 'hydropathy': -1.0,
          'aliphaticity': -0.157024793388, 'hydrogenation': -0.401797175866, 'aromaticity': -0.0642673521851,
          'hydroxyethilation': -0.514863258026}
}
