from unittest import TestCase, skipUnless
from six import assertRaisesRegex

from dark.aa import CODONS
from dark.diamond.run import DiamondExecutor, diamondInstalled
from dark.reads import Read, Reads

from .sample_proteins import SAMPLE_DATA


@skipUnless(diamondInstalled(), 'DIAMOND is not installed')
class TestDiamondExecutor(TestCase):
    """
    Test the execution of DIAMOND.
    """
    def testEmptyDatabase(self):
        """
        If no subject sequences have been added, doing a search must result in
        a ValueError.
        """
        de = DiamondExecutor()
        queries = Reads([Read('id', 'ACGT')])
        error = '^No subject sequences in the database$'
        assertRaisesRegex(self, ValueError, error, list, de.search(queries))
        de.cleanup()

    def testNoQueries(self):
        """
        If no query sequences are passed, doing a search must result in
        a ValueError.
        """
        de = DiamondExecutor()
        de.addSubject(Read('id', 'ACGT'))
        queries = Reads()
        error = '^No query sequences were passed$'
        assertRaisesRegex(self, ValueError, error, list, de.search(queries))
        de.cleanup()

    def testOne(self):
        """
        If one query sequence is passed, doing a search must result in the
        expected result.
        """
        qid = 'query'
        qseq = 'CCAGAGCGCACATACTGAATAGCAAACGGATTCTCCTTCGGGTCACACTCAATTGGG'
        qqual = 'EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE'
        sid = 'subject'

        with DiamondExecutor() as de:
            de.addSubject(Read(
                sid,
                'MVMWMEVPCSGVDKYIDSLEKFIEKRIPVKMTGKRNNKLRVKIVAKKATKKKNEVTKLG'
                'AALRALGGLGGGAVGSLFGNPATGSTIGTGLGAALSRWLGRGDYRVSSNSVVQQSLKGT'
                'SSIPSMHQDGQSVVVRHKEFVTEVRGNTSFLVRGTFDINPGRAETFPWLAGVASRFQEY'
                'KIRGLVWHYVPSSGTAVSGTDPALGTVMLQTSYRSNDIPPANKMEVLNEYWSSESVPSE'
                'AFCHPIECDPKENPFNIQYVRTDKVPDGDSKLLYDLGTTHLCVSGQQANDVVLGDLWCT'
                'YEIELKKPIVNSNVTSVARSAALAFTGTLDLNSWFNGSVVNFGSLDVTANVKTISFPAR'
                'LTGRFLVTVTIIASTTFTAADLSGTPLTTNCSLFALPTGSTYLRTVLGGTTPTVNIATL'
                'QFAVDIQDSAKTASVTCLGSFTGAATRSMVTVTPYLM'))
            queries = Reads([Read(qid, qseq, qqual)])
            (result,) = list(de.search(queries))

        self.assertEqual(
            {
                'bitscore': 40.0,
                'btop': '11AN5ST',
                'qframe': -2,
                'qend': 3,
                'full_qqual': qqual,
                'qlen': 57,
                'full_qseq': qseq,
                'qseqid': qid,
                'qstart': 56,
                'slen': 450,
                'sstart': 241,
                'stitle': sid,
            },
            result)

    def testYP_009259545(self):
        """
        Test for a match against YP_009259545
        """
        proteinAccession = 'YP_009259545.1'
        proteinSequence = SAMPLE_DATA['proteins'][proteinAccession]['protein']
        proteinId = SAMPLE_DATA['proteins'][proteinAccession]['id']

        qid = 'query'
        qseq = ''.join(CODONS[aa][0] for aa in proteinSequence[10:50])
        qqual = 'E' * len(qseq)

        with DiamondExecutor() as de:
            de.addSubject(Read(proteinId, proteinSequence))
            queries = Reads([Read(qid, qseq, qqual)])
            (result,) = list(de.search(queries))

        self.assertEqual(
            {
                'bitscore': 82.4,
                'btop': '40',
                'qframe': 1,
                'qend': 120,
                'full_qqual': qqual,
                'qlen': len(qseq),
                'full_qseq': qseq,
                'qseqid': 'query',
                'qstart': 1,
                'slen': len(proteinSequence),
                'sstart': 11,
                'stitle': proteinId,
            },
            result)
