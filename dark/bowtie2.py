import multiprocessing
from tempfile import mkdtemp
import os
from os.path import join
import requests

from dark.process import Executor


class Bowtie2(object):
    """
    Run Bowtie2.
    """

    def __init__(self, executor=None, threads=None, verboseFp=None,
                 dryRun=False, reference=None, tempdir=None,
                 tmpChmod=None):
        self._executor = executor or Executor(dryRun)
        if dryRun:
            self.tempdir = tempdir or '/tmp/xxx'
        else:
            self.tempdir = mkdtemp(prefix='bt2-', dir=tempdir)
            if tmpChmod:
                self._executor.execute(f'chmod {tmpChmod} {self.tempdir}')
        self._samFile = join(self.tempdir, 'result.sam')
        self._bamFile = join(self.tempdir, 'result.bam')
        self._indexFile = join(self.tempdir, 'index')
        self._verboseFp = verboseFp
        self._indexCalled = self._samExists = self._bamExists = False
        self.nThreads = threads or multiprocessing.cpu_count()
        self._reference = reference

    def buildIndex(self, index):
        """
        Find or make a Bowtie2 index.
        """
        if os.path.exists(index):
            # Check if this is a pre-existing bowtie2 index. Look for and
            # remove a bowtie2 suffix, if any (otherwise bowtie2 will
            # complain). We do things this way to allow a user to use TAB
            # completion on the command line to give them the full path to
            # any bowtie index file.
            for suffix in ('1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 '
                           'rev.2.bt2').split():
                suffix = '.' + suffix
                if index.endswith(suffix):
                    self._indexFile = index[:-len(suffix)]
                    self._report('Using pre-existing Bowtie2 index %r.' %
                                 self._indexFile)
                    break
            else:
                # Assume a FASTA file and make an index.
                self._indexFile = self._makeIndexFromFastaFile(index)
        else:
            # Not a filename. So either the start of the path to a bowtie2
            # index or else an accession number.
            if os.path.exists(index + '.1.bt2'):
                self._report('Using pre-existing Bowtie2 index %r.' % index)
                self._indexFile = index
            else:
                # Assume an accession number.
                self._indexFile = self._makeIndexFromAccession(index)

        self._indexCalled = True

    def align(self, bowtie2Args='--no-unal', fastq1=None, fastq2=None,
              threads=None, discardSAM=False, readGroup='orig',
              sampleName='orig'):
        """
        Run Bowtie2.
        """
        if not self._indexCalled:
            raise ValueError('buildIndex() has not yet been called.')

        self._report('Aligning with Bowtie2.')

        nThreads = threads or self.nThreads
        samFile = '/dev/null' if discardSAM else self._samFile

        if fastq1 and fastq2:
            self._executor.execute(
                "bowtie2 %s --threads %d --rg-id '%s' --rg 'SM:%s' -x '%s' "
                "-1 '%s' -2 '%s' > '%s'" % (
                    bowtie2Args, nThreads, readGroup, sampleName,
                    self._indexFile, fastq1, fastq2, samFile))
        elif fastq1:
            self._executor.execute(
                "bowtie2 %s --threads %d --rg-id '%s' --rg 'SM:%s' -x '%s' "
                "-U '%s' > '%s'" % (
                    bowtie2Args, nThreads, readGroup, sampleName,
                    self._indexFile, fastq1, samFile))
        else:
            raise ValueError('At least fastq1 must be passed.')

        self._samExists = True

    def _report(self, mesg):
        """
        Print a progress message.
        """
        if self._verboseFp:
            print(mesg, file=self._verboseFp)

    def makeBAM(self, samtoolsViewArgs=''):
        """
        Convert the SAM to BAM.
        """
        # The following will raise if there is no SAM file.
        self._SAMorBAM()
        self._report('Converting SAM to BAM.')
        self._executor.execute(
            "samtools view -b %s '%s' > '%s'" %
            (samtoolsViewArgs, self._samFile, self._bamFile))
        self._bamExists = True

    def _SAMorBAM(self):
        """
        Do we have a SAM or BAM file available?
        """
        if self._bamExists:
            return 'BAM'
        elif self._samExists:
            return 'SAM'
        else:
            raise ValueError('bowtie2() has not yet been called.')

    def outputFile(self):
        """
        The name of the output file.
        """
        if self._bamExists:
            return self._bamFile
        elif self._samExists:
            return self._samFile
        else:
            raise ValueError('bowtie2() has not yet been called.')

    def indexBAM(self):
        """
        Index the BAM file.
        """
        which = self._SAMorBAM()

        if which != 'BAM':
            raise ValueError('makeBAM() has not yet been called.')

        self._report('Indexing BAM.')
        self._executor.execute("samtools index '%s'" % self._bamFile)

    def sort(self, byName=False):
        """
        Sort the BAM or SAM.
        """
        which = self._SAMorBAM()
        self._report('Sorting %s (by %s).' % (
            which, 'name' if byName else 'coord'))
        inFile = self._bamFile if which == 'BAM' else self._samFile
        sortedFile = join(self.tempdir, 'result-sorted.' + which.lower())
        self._executor.execute(
            "samtools sort %s'%s' > '%s'" % (
                '-n ' if byName else '', inFile, sortedFile))
        self._executor.execute("mv '%s' '%s'" % (sortedFile, inFile))

    def removePrimers(self, bedFile):
        """
        Removes primers specified in the bed file
        """
        which = self._SAMorBAM()

        if which != 'BAM':
            raise ValueError('makeBAM() has not yet been called.')

        self._report("removing primers specified in %s" % bedFile)
        tempTrimmedBamPrefix = "%s.trimmed" % self._bamFile
        self._executor.execute(
            "ivar trim -b '%s' -p '%s' -i '%s' -q 20 -m 30 -s 4 -e" %
            (bedFile, tempTrimmedBamPrefix, self._bamFile))
        self._executor.execute("mv '%s'.bam '%s'" %
                               (tempTrimmedBamPrefix, self._bamFile))

    def markDuplicatesPicard(self, picardFile):
        """
        Use Picard to mark duplicates.
        """
        which = self._SAMorBAM()
        self._report('Marking duplicates with Picard.')

        inFile = self._bamFile if which == 'BAM' else self._samFile
        tempFile = join(self.tempdir, 'picard-duplicates.' + which.lower())
        tempErrFile = join(self.tempdir, 'picard.errs')

        self._executor.execute(
            'java -Xmn2g -Xms2g -Xmx2g -jar %s '
            "MarkDuplicates I='%s' O='%s' M=/dev/null >'%s' 2>&1" %
            (picardFile, inFile, tempFile, tempErrFile))

        self._executor.execute("mv '%s' '%s'" % (tempFile, inFile))

    def markDuplicatesGATK(self, threads=None):
        """
        Use GATK to mark duplicates.
        """
        nThreads = threads or self.nThreads
        which = self._SAMorBAM()
        self._report('Marking duplicates with GATK.')

        inFile = self._bamFile if which == 'BAM' else self._samFile
        tempFile = join(self.tempdir, 'gatk-duplicates.' + which.lower())

        self._executor.execute(
            'gatk MarkDuplicatesSpark '
            "-I '%s' -O '%s' --conf spark.executor.cores=%d" %
            (inFile, tempFile, nThreads))

        self._executor.execute("mv '%s' '%s'" % (tempFile, inFile))

    def callHaplotypesGATK(self, picardJar, vcfFile=None, referenceFasta=None):
        """
        Use GATK to call haplotypes.
        """
        which = self._SAMorBAM()
        self._report('Calling haplotypes with GATK.')

        inFile = self._bamFile if which == 'BAM' else self._samFile
        vcfFile = vcfFile or join(self.tempdir, 'output.vcf.gz')

        if referenceFasta is None:
            if self._reference:
                self._report('Using %s as a reference.' % self._reference)
                referenceFasta = self._reference
            else:
                raise ValueError('No reference was passed, given to the '
                                 'Bowtie2 __init__, or used in buildIndex.')

        indexFile = referenceFasta + '.fai'
        if os.path.exists(indexFile):
            removeIndex = False
        else:
            removeIndex = True
            self._executor.execute("samtools faidx '%s'" % referenceFasta)

        if referenceFasta.lower().endswith('.fasta'):
            dictFile = referenceFasta[:-len('.fasta')] + '.dict'
        else:
            dictFile = referenceFasta + '.dict'

        if os.path.exists(dictFile):
            removeDict = False
        else:
            removeDict = True
            self._executor.execute(
                "java -jar '%s' CreateSequenceDictionary R='%s' O='%s'" %
                (picardJar, referenceFasta, dictFile))

        self._executor.execute(
            'gatk --java-options -Xmx4g HaplotypeCaller '
            "--reference '%s' "
            "--input '%s' "
            "--output '%s' "
            "--sample-ploidy 1 "
            "--dont-use-soft-clipped-bases true "
            '-ERC GVCF' %
            (referenceFasta, inFile, vcfFile))

        if removeIndex:
            self._executor.execute("rm '%s'" % indexFile)

        if removeDict:
            self._executor.execute("rm '%s'" % dictFile)

    def callHaplotypesBcftools(self, vcfFile=None, referenceFasta=None):
        """
        Use bcftools call to call haplotypes.
        """
        which = self._SAMorBAM()
        self._report('Calling haplotypes with bcftools call.')

        inFile = self._bamFile if which == 'BAM' else self._samFile
        vcfFile = vcfFile or join(self.tempdir, 'output.vcf.gz')

        if referenceFasta is None:
            if self._reference:
                self._report('Using %s as a reference.' % self._reference)
                referenceFasta = self._reference
            else:
                raise ValueError('No reference was passed, given to the '
                                 'Bowtie2 __init__, or used in buildIndex.')

        self._executor.execute(
            'bcftools mpileup --max-depth 5000 -Ou -f "%s" "%s" | '
            'bcftools call --ploidy 1 -mv -Oz -o "%s"' %
            (referenceFasta, inFile, vcfFile))

        self._executor.execute('bcftools index %s' % vcfFile)

    def removeDuplicates(self):
        """
        Use samtools to remove marked duplicates.
        """
        which = self._SAMorBAM()
        self._report('Removing marked duplicates.')

        inFile = self._bamFile if which == 'BAM' else self._samFile
        tempFile = join(self.tempdir, 'non-duplicates.' + which.lower())

        self._executor.execute(
            "samtools view -b -F 1024 '%s' > '%s'" % (inFile, tempFile))
        self._executor.execute("mv '%s' '%s'" % (tempFile, inFile))

    def _makeIndexFromFastaFile(self, fastaFilename):
        index = join(self.tempdir, 'index')

        self._report('Building Bowtie2 index from %s.' % fastaFilename)
        self._executor.execute("bowtie2-build --quiet '%s' '%s'" %
                               (fastaFilename, index))

        if self._reference is None:
            self._report('Will use this FASTA file as the reference '
                         '(if needed).')
            self._reference = fastaFilename

        return index

    def _makeIndexFromAccession(self, accessionId):
        fastaFilename = join(self.tempdir, accessionId + '.fasta')
        index = join(self.tempdir, 'index')

        self._report('Downloading FASTA for accession %s from NCBI.' %
                     accessionId)

        URL = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
               'db=nucleotide&id=%s&rettype=fasta&retmode=text' % accessionId)

        if not self._executor.dryRun:
            with open(fastaFilename, 'w') as fp:
                print(requests.get(URL).text.rstrip('\n'), file=fp)

        self._report('Building bowtie2 index for %s.' % accessionId)
        self._executor.execute("bowtie2-build --quiet '%s' '%s'" %
                               (fastaFilename, index))

        return index

    def close(self):
        """
        Clean up.
        """
        self._executor.execute("rm -r '%s'" % self.tempdir)
