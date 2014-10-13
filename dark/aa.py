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
