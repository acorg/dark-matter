import subprocess
from cStringIO import StringIO
from Bio import SeqIO


def getSequence(hitId, db='nt'):
    """
    hitId: the hit_id field from a BLAST record hsp. Of the form
        'gi|63148399|gb|DQ011818.1|' or anything recognized by the -entry param
        of blastdbcmd.
    db: the str name of the BLAST database to search.
    """
    fasta = subprocess.check_output(
        ['blastdbcmd', '-entry', hitId, '-db', db])
    return SeqIO.read(StringIO(fasta), 'fasta')
