from Bio import Entrez, SeqIO
from six.moves.urllib_error import URLError

Entrez.email = 'tcj25@cam.ac.uk'


def getSequence(title, db='nucleotide'):
    """
    Get information about a sequence from Genbank.

    @param title: A C{str} sequence title from a BLAST hit. Of the form
        'gi|63148399|gb|DQ011818.1| Description...'.
    @param db: The C{str} name of the Entrez database to consult.

    NOTE: this uses the network!  Also, there is a 3 requests/second limit
    imposed by NCBI on these requests so be careful or your IP will be banned.
    """
    titleId = title.split(' ', 1)[0]
    try:
        gi = titleId.split('|')[1]
    except IndexError:
        # Assume we have a gi number directly, and make sure it's a string.
        gi = str(titleId)

    try:
        client = Entrez.efetch(db=db, rettype='gb', retmode='text', id=gi)
    except URLError:
        return None
    else:
        record = SeqIO.read(client, 'gb')
        client.close()
        return record
