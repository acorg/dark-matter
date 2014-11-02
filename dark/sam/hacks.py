from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess as sp


def checkSAMfile(samFilename):
    """
    Checks that the file inputted as a SAM file is indeed a SAM file.
    NB header section is optional but length of subjects is stored
    there - so it is compulsory.
    @param samFilename: a C{str} of a SAM file.
    """
    headerLines = 0
    nonHeaderLines = 0
    with open(samFilename) as samFile:
        for line in samFile:
            while line[0] == '@':
                headerLines += 1
                break
            else:
                assert line[0] != '@', 'A header line is in alignment section'
                elements = line.strip().split()
                assert len(elements) > 10, ('SAM file %s does not contain '
                                            'at least 11 fields.'
                                            % samFilename)
                nonHeaderLines += 1
    assert headerLines > 0 or nonHeaderLines > 0, ('SAM file %s is empty.'
                                                   % samFilename)
    assert headerLines > 0, ('SAM file %s does not contain header section '
                             'which is required.' % samFilename)
    assert nonHeaderLines > 0, ('SAM file %s does not contain alignment '
                                'section.' % samFilename)
    return True


def checkFASTAfile(fastaFilename):
    """
    Checks that the file inputted as a FASTA file is indeed a FASTA file.
    @param fastaFilename: a C{str} of a FASTA file.
    """
    with open(fastaFilename) as fastaFile:
        # Make into an iterable so can compare two lines.
        fastaFile = iter(fastaFile)
        for line in fastaFile:
            assert line[0] == '>', 'FASTA file does not begin with a title'
            line = next(fastaFile)
            # Invalid if two lines begin with >
            assert line[0] != '>' and line, 'Invalid FASTA file format'


def checkFASTQfile(fastqFilename):
    """
    TODO: Finish
    Checks that the file inputted as a FASTQ file is indeed a FASTQ file.
    @param fastqFilename: a C{str} of a FASTQ file.
    """
    with open(fastqFilename) as fastqFile:
        fastqFile = iter(fastqFile)
        for line in fastqFile:
            header = line[1:]
            assert line[0] == '@', 'Invalid header of entry: %s' % header
            line = next(fastqFile)
            assert line, 'Empty raw sequence in entry: %s' % header
            line = next(fastqFile)
            assert line[0] == '+', 'Invalid third line of entry: %s' % header
            line = next(fastqFile)
            assert line, 'Empty quality score in entry: %s' % header
    return True


def samSubtract(samFile, outFile):
    """
    Takes a SAM file, makes a set of the seqids of unaligned sequences.
    Makes a new FASTA or FASTQ file with just these seqs (outFile).

    @param samFile: a C{str} name of a SAM file.
    @param outFile: a C{str} name of the output file which should be a .fasta
        or .fastq file.
    """
    if not checkFASTAfile(outFile) and not outFile.lower().endswith('.fastq'):
        raise ValueError('An output FASTA or FASTQ file must be given')
    if checkSAMfile(samFile):
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
                ofp.write('@%s\n%s\n+\n%s\n' % (item['query'], item['seq'],
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


def findMD(samFile, fastaFile):
    """
    Calculates the optional MD tag. If MD tag already present
    raises error if any are found to be different.
    Requires samtools.

    @param samFile: a C{str} of a SAM file.
    @param fastaFile: a C{str} of the FASTA file used for the
        SAM file.
    @return: samFile with the MD strings.
    """
    if checkSAMfile(samFile) and checkFASTAfile(fastaFile):
        samFileNew = ''.join(samFile.split())[:-4] + '-fillmd.sam'
        sp.Popen(['samtools', 'fillmd', '-S', samFile, fastaFile, '>',
                 samFileNew], stderr=sp.PIPE)
        return samFileNew
