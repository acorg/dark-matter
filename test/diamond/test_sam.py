from unittest import TestCase
from six import StringIO

from dark.aa import CODONS
from dark.diamond.sam import SimpleDiamondSAMWriter
from dark.diamond.run import DiamondExecutor
from dark.proteins import SqliteIndex
from dark.reads import Read, Reads

from .sample_proteins import SAMPLE_DATA


class TestSimpleDiamondSAMWriter(TestCase):

    def testTibetanFrogHBV(self):

        proteinAccession = 'YP_009259545.1'
        proteinSequence = SAMPLE_DATA['proteins'][proteinAccession]['protein']
        proteinId = SAMPLE_DATA['proteins'][proteinAccession]['id']
        proteinRange = SAMPLE_DATA['proteins'][proteinAccession]['range']
        queryStartInProtein = 10  # This is a 0-based amino acid offset.
        queryLenInProtein = 40  # This is an amino acid length.

        genomeAccession = 'NC_030446.1'
        genomeSequence = SAMPLE_DATA['genomes'][genomeAccession]['genome']

        # The query sequence is nucleotides that match the amino acids in the
        # protein. Here we use the first ([0] in the below) codon for each
        # amino acid to make the nucleotide sequence.
        queryId = 'query'
        querySequence = ''.join(
            CODONS[aa][0]
            for aa in proteinSequence[queryStartInProtein:
                                      queryStartInProtein + queryLenInProtein])
        queryQuality = 'E' * len(querySequence)

        # Use the protein sequence to make a DIAMOND database and run DIAMOND
        # on the query. Yes, this really runs DIAMOND, so you need to have it
        # installed, with its executable somewhere in your shell's PATH.
        with DiamondExecutor() as de:
            de.addSubject(Read(proteinId, proteinSequence))
            queries = Reads([Read(queryId, querySequence, queryQuality)])
            (diamondResult,) = list(de.search(queries))

        # Make sure DIAMOND gives us back what we expected.
        self.assertEqual(
            {
                'bitscore': 83.6,
                'btop': str(queryLenInProtein),  # Exact match of all AAs.
                'qframe': 1,
                'qend': 3 * queryLenInProtein,
                'full_qqual': queryQuality,
                'qlen': len(querySequence),
                'full_qseq': querySequence,
                'qseqid': 'query',
                'qstart': 1,
                'slen': len(proteinSequence),
                'sstart': queryStartInProtein + 1,  # DIAMOND is 1-based.
                'stitle': proteinId,
            },
            diamondResult)

        # Make a genomes/proteins sqlite database and add information about
        # the protein and the nucleotide genome it comes from.
        db = SqliteIndex(':memory:')
        with db._connection as connection:
            connection.execute(
                'INSERT INTO '
                'genomes(accession, sequence, proteinCount) '
                'VALUES (?, ?, ?)',
                (genomeAccession, genomeSequence, 5))

            connection.execute(
                'INSERT INTO '
                'proteins(accession, sequence, offsets, reversed, circular, '
                'gene, note) VALUES (?, ?, ?, ?, ?, ?, ?)',
                (proteinAccession, proteinSequence, proteinRange, 0, 0, None,
                 None))

        # Make a DIAMOND to SAM writer and give it the DIAMOND output.
        writer = SimpleDiamondSAMWriter(db)
        writer.addMatch(
            '\t'.join(map(str, (
                diamondResult['bitscore'],
                diamondResult['btop'],
                diamondResult['qframe'],
                diamondResult['qend'],
                diamondResult['full_qqual'],
                diamondResult['qlen'],
                diamondResult['full_qseq'],
                diamondResult['qseqid'],
                diamondResult['qstart'],
                diamondResult['slen'],
                diamondResult['sstart'],
                diamondResult['stitle'],
            ))))

        # Tell the writer to save the matches as SAM and check the result.
        fp = StringIO()
        writer.save(filename=fp)

        flags = '31'
        self.assertEqual(
            '\n'.join([
                '@SQ\tSN:%s\tLN:%d' % (genomeAccession, len(genomeSequence)),
                '\t'.join([
                    queryId,
                    '0',
                    genomeAccession,
                    flags,
                    '255',
                    '120M',  # (Exact) match of 40 AAs.
                    '*',
                    '0',
                    '0',
                    querySequence,
                    queryQuality,
                    'AS:i:%d' % int(diamondResult['bitscore']),
                ])
            ]) + '\n',
            fp.getvalue())
