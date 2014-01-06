from Bio import Entrez, SeqIO
from urllib2 import URLError

Entrez.email = 'tcj25@cam.ac.uk'


def getSequence(hitId, db='nucleotide'):
    """
    Get information about a sequence from Genbank.

    hitId: either a hit from a BLAST record, in the form
        'gi|63148399|gb|DQ011818.1|' in which case we use the 2nd field, as
        delimited by '|', to fetch from Genbank.  Or, a gi number (the 2nd
        field just mentioned).
    @param db: The C{str} name of the Genbank database to consult.

    NOTE: this uses the network!  Also, there is a 3 requests/second limit
    imposed by NCBI on these requests so be careful or your IP will be banned.
    """
    try:
        gi = hitId.split('|')[1]
    except IndexError:
        # Assume we have a gi number directly, and make sure it's a string.
        gi = str(hitId)

    try:
        client = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text',
                               id=gi)
    except URLError:
        return None
    else:
        record = SeqIO.read(client, 'gb')
        client.close()
        return record
