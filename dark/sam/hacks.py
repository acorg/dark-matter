from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess as sp


def samSubtract(samFile, outFile):
    """
    Takes a SAM file, makes a set of the seqids of unaligned sequences.
    Makes a new FASTA or FASTQ file with just these seqs (outFile).

    @param samFile: a C{str} name of a SAM file.
    @param outFile: a C{str} name of the output file which should be a .fasta
        or .fastq file.
    """
    # Put these in bin script once written
    if not samFile.lower().endswith('.sam'):
        raise ValueError('A SAM file must be given')
    if (not outFile.lower().endswith('.fasta') and not
            outFile.lower().endswith('.fastq')):
        raise ValueError('An output FASTA or FASTQ file must be given')

    records = []
    with open(samFile) as fp:
        for line in fp:
            if not line.startswith('$') and not line.startswith('@'):
                line = line.strip().split()
                refSeqName = line[2]
                if refSeqName == '*':
                    query = line[0]
                    seq = line[9]
                    if outFile.lower().endswith('.fastq'):
                        qual = line[10]
                        record = {'query': query, 'seq': seq, 'qual': qual}
                        # record = SeqRecord(Seq(seq), id=query)
                        # record.letter_annotations['phred_quality']=qual
                        records.append(record)
                    else:
                        record = SeqRecord(Seq(seq), id=query)
                        records.append(record)

    with open(outFile, 'w') as ofp:
        if outFile.lower().endswith('.fasta'):
            # NB seq in fasta may be over multiple lines
            SeqIO.write(records, ofp, 'fasta')
        else:
            # Previously: SeqIO.write(records, ofp, 'fastq')
            # The SeqIO.write for fastq requires you to define the type of
            # quality string eg illumina, solexa etc. See help on QualityIO.
            # https://github.com/biopython/biopython/blob/master/Bio/SeqIO/QualityIO.py
            # Need to use FastqGeneralIterator? But that req file input.

            for item in records:
                ofp.write("@%s\n%s\n+\n%s\n" % (item['query'], item['seq'],
                                                item['qual']))


def bamSubtract(bamFile, fastaFile, outFile):
    """
    Takes a BAM file, converts to a SAM file using samtools and then calls
    the samSubtract function. Useful if BAM files are stored instead of SAM
    files. Samtools required.

    @param bamFile: a C{str} name of a BAM file.
    @param fastaFile: a C{str} name of a FASTA file or FASTQ file.
    @param outFile: a C{str} name of the output file which
        should be a .fasta or .fastq file.
    """
    # Need to put in bin file
    if not bamFile.lower().endswith('.bam'):
        raise ValueError('A BAM file must be given')
    if (not fastaFile.lower().endswith('.fasta') and not
            fastaFile.lower().endswith('.fastq')):
        raise ValueError('A FASTA or FASTQ file must be given')
    if (not outFile.lower().endswith('.fasta') and not
            outFile.lower().endswith('.fastq')):
        raise ValueError('An output FASTA or FASTQ file must be given')

    samFileNew = ''.join(bamFile.split())[:-3] + 'sam'
    sp.Popen(['samtools', 'view', '-h', bamFile, '>', samFileNew],
             stderr=sp.PIPE)

    return samSubtract(samFileNew, fastaFile, outFile)
