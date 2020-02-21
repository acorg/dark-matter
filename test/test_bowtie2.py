from unittest import TestCase
from io import StringIO

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from dark.bowtie2 import Bowtie2
from dark.process import Executor


class TestBowtie2(TestCase):
    """
    Test the Bowtie2 class.
    """
    def testIndexAccession(self):
        """
        Making a Bowtie index from an accession number must result in the
        expected commands being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('MN908947.3')
        self.assertEqual(
            "$ bowtie2-build --quiet '/tmp/xxx/MN908947.3.fasta' "
            "'/tmp/xxx/index'",
            e.log[-1])
        log = fp.getvalue()
        self.assertTrue(
            log.startswith('Downloading FASTA for accession MN908947.3 '
                           'from NCBI.\n'))

    @patch('os.path.exists')
    def testIndexFromFile(self, existsMock):
        """
        Making a Bowtie index from a file name must result in the
        expected commands being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        self.assertEqual(
            "$ bowtie2-build --quiet 'file.fasta' '/tmp/xxx/index'",
            e.log[-1])
        log = fp.getvalue()
        self.assertTrue(
            log.find("Building Bowtie2 index from file.fasta.\n") > -1)

    @patch('os.path.exists')
    def testIndexFromBowtie2File(self, existsMock):
        """
        Making a Bowtie index from a file that is an existing Bowtie2 index
        file must result in the expected commands being run and output being
        produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.1.bt2')
        self.assertEqual(-1, e.log[-1].find('bowtie2-build'))
        log = fp.getvalue()
        self.assertEqual("Using pre-existing Bowtie2 index 'file'.\n", log)

    @patch('os.path.exists')
    def testIndexFromBowtie2FilePrefix(self, existsMock):
        """
        Making a Bowtie index from a file that is a prefix of an existing
        Bowtie2 index file name must result in the expected commands being run
        and output being produced.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename):
                if self.count == 0:
                    self.test.assertEqual('idx-file', filename)
                    self.count += 1
                    return False
                elif self.count == 1:
                    self.test.assertEqual('idx-file.1.bt2', filename)
                    self.count += 1
                    return True
                else:
                    self.test.fail('Unexpected third call to exists.')

        existsMock.side_effect = SideEffect(self).sideEffect

        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('idx-file')
        self.assertEqual(-1, e.log[-1].find('bowtie2-build'))
        log = fp.getvalue()
        self.assertEqual("Using pre-existing Bowtie2 index 'idx-file'.\n", log)

    @patch('os.path.exists')
    def testAlignOneFASTQ(self, existsMock):
        """
        Making a Bowtie index from a file name and running an alignment with
        one FASTQ file must result in the expected commands being run and
        output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq', threads=4)
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nAligning with Bowtie2.\n'))
        self.assertEqual(
            "$ bowtie2 --no-unal --threads 4 --rg-id 'orig' --rg 'SM:orig' "
            "-x '/tmp/xxx/index' "
            "-U 'file1.fastq' > '/tmp/xxx/result.sam'", e.log[-1])

    @patch('os.path.exists')
    def testAlignOneFASTQWithReadGroup(self, existsMock):
        """
        Making a Bowtie index from a file name and running an alignment with
        one FASTQ file and specifying a read group must result in the expected
        commands being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq', threads=4, readGroup='xxx')
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nAligning with Bowtie2.\n'))
        self.assertEqual(
            "$ bowtie2 --no-unal --threads 4 --rg-id 'xxx' --rg 'SM:orig' "
            "-x '/tmp/xxx/index' "
            "-U 'file1.fastq' > '/tmp/xxx/result.sam'", e.log[-1])

    @patch('os.path.exists')
    def testAlignOneFASTQWithSampleName(self, existsMock):
        """
        Making a Bowtie index from a file name and running an alignment with
        one FASTQ file and specifying a sample name must result in the expected
        commands being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq', threads=4, sampleName='xxx')
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nAligning with Bowtie2.\n'))
        self.assertEqual(
            "$ bowtie2 --no-unal --threads 4 --rg-id 'orig' --rg 'SM:xxx' "
            "-x '/tmp/xxx/index' "
            "-U 'file1.fastq' > '/tmp/xxx/result.sam'", e.log[-1])

    @patch('os.path.exists')
    def testAlignOneFASTQWithFlags(self, existsMock):
        """
        Making a Bowtie index from a file name and running an alignment with
        one FASTQ file and some (fictional) args must result in the expected
        commands being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq', threads=4, bowtie2Args='--up --dn')
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nAligning with Bowtie2.\n'))
        self.assertEqual(
            "$ bowtie2 --up --dn --threads 4 --rg-id 'orig' --rg 'SM:orig' -x "
            "'/tmp/xxx/index' -U 'file1.fastq' > '/tmp/xxx/result.sam'",
            e.log[-1])

    @patch('os.path.exists')
    def testAlignTwoFASTQs(self, existsMock):
        """
        Making a Bowtie index from a file name and running an alignment with
        two FASTQ files must result in the expected commands being run and
        output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq', fastq2='file2.fastq', threads=4)
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nAligning with Bowtie2.\n'))
        self.assertEqual(
            "$ bowtie2 --no-unal --threads 4 --rg-id 'orig' --rg 'SM:orig' -x "
            "'/tmp/xxx/index' "
            "-1 'file1.fastq' -2 'file2.fastq' > '/tmp/xxx/result.sam'",
            e.log[-1])

    @patch('os.path.exists')
    def testAlignTwoFASTQsWithReadGroup(self, existsMock):
        """
        Making a Bowtie index from a file name and running an alignment with
        two FASTQ files and a read group must result in the expected commands
        being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq', fastq2='file2.fastq', threads=4,
                 readGroup='xxx')
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nAligning with Bowtie2.\n'))
        self.assertEqual(
            "$ bowtie2 --no-unal --threads 4 --rg-id 'xxx' --rg 'SM:orig' -x "
            "'/tmp/xxx/index' "
            "-1 'file1.fastq' -2 'file2.fastq' > '/tmp/xxx/result.sam'",
            e.log[-1])

    @patch('os.path.exists')
    def testAlignTwoFASTQsWithSampleName(self, existsMock):
        """
        Making a Bowtie index from a file name and running an alignment with
        two FASTQ files and a sample name must result in the expected commands
        being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq', fastq2='file2.fastq', threads=4,
                 sampleName='xxx')
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nAligning with Bowtie2.\n'))
        self.assertEqual(
            "$ bowtie2 --no-unal --threads 4 --rg-id 'orig' --rg 'SM:xxx' -x "
            "'/tmp/xxx/index' "
            "-1 'file1.fastq' -2 'file2.fastq' > '/tmp/xxx/result.sam'",
            e.log[-1])

    @patch('os.path.exists')
    def testMakeBAM(self, existsMock):
        """
        Making a BAM file must result in the expected commands being run and
        output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq')
        bt.makeBAM()
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nConverting SAM to BAM.\n'))
        self.assertEqual(
            "$ samtools view -b  '/tmp/xxx/result.sam' > "
            "'/tmp/xxx/result.bam'",
            e.log[-1])

    @patch('os.path.exists')
    def testMakeIndexedBAM(self, existsMock):
        """
        Making a BAM file and indexing it must result in the expected commands
        being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq')
        bt.makeBAM()
        bt.indexBAM()
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nIndexing BAM.\n'))
        self.assertEqual("$ samtools index '/tmp/xxx/result.bam'", e.log[-1])

    @patch('os.path.exists')
    def testRemoveDuplicates(self, existsMock):
        """
        Removing duplicates must result in the expected commands
        being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq')
        bt.removeDuplicates()
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nRemoving marked duplicates.\n'))
        self.assertEqual("$ samtools view -b -F 1024 '/tmp/xxx/result.sam' "
                         "> '/tmp/xxx/non-duplicates.sam'", e.log[-2])
        self.assertEqual("$ mv '/tmp/xxx/non-duplicates.sam' "
                         "'/tmp/xxx/result.sam'", e.log[-1])

    @patch('os.path.exists')
    def testSortByCoord(self, existsMock):
        """
        Sorting by coord (the default) must result in the expected commands
        being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq')
        bt.sort()
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nSorting SAM (by coord).\n'))
        self.assertEqual("$ samtools sort '/tmp/xxx/result.sam' > "
                         "'/tmp/xxx/result-sorted.sam'", e.log[-2])
        self.assertEqual("$ mv '/tmp/xxx/result-sorted.sam' "
                         "'/tmp/xxx/result.sam'", e.log[-1])

    @patch('os.path.exists')
    def testSortByName(self, existsMock):
        """
        Sorting by name must result in the expected commands
        being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq')
        bt.sort(byName=True)
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nSorting SAM (by name).\n'))
        self.assertEqual("$ samtools sort -n '/tmp/xxx/result.sam' > "
                         "'/tmp/xxx/result-sorted.sam'", e.log[-2])
        self.assertEqual("$ mv '/tmp/xxx/result-sorted.sam' "
                         "'/tmp/xxx/result.sam'", e.log[-1])

    @patch('os.path.exists')
    def testOutputFileSAM(self, existsMock):
        """
        The output file must be a SAM file if BAM hasn't been produced.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename):
                if self.count == 0:
                    self.test.assertEqual('index', filename)
                    self.count += 1
                    return False
                elif self.count == 1:
                    self.test.assertEqual('index.1.bt2', filename)
                    self.count += 1
                    return True
                elif self.count == 2:
                    self.test.assertEqual('/tmp/xxx/result.bam', filename)
                    self.count += 1
                    return False
                elif self.count == 3:
                    self.test.assertEqual('/tmp/xxx/result.sam', filename)
                    self.count += 1
                    return True
                else:
                    self.test.fail(
                        'Unexpected 5th call to exists. Filename: %r.' %
                        filename)

        existsMock.side_effect = SideEffect(self).sideEffect

        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('index')
        bt.align(fastq1='file1.fastq')
        self.assertEqual('/tmp/xxx/result.sam', bt.outputFile())

    @patch('os.path.exists')
    def testOutputFileBAM(self, existsMock):
        """
        The output file must be BAM if it has been produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq')
        bt.makeBAM()
        self.assertEqual('/tmp/xxx/result.bam', bt.outputFile())

    @patch('os.path.exists')
    def testMarkDuplicatesPicard(self, existsMock):
        """
        Using Picard to mark duplicates must result in the expected commands
        being run and output being produced.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq')
        bt.makeBAM()
        bt.indexBAM()
        bt.markDuplicatesPicard('picard.jar')
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nMarking duplicates with Picard.\n'))
        self.assertEqual(
            '$ java -Xmn2g -Xms2g -Xmx2g -jar picard.jar MarkDuplicates '
            "I='/tmp/xxx/result.bam' O='/tmp/xxx/picard-duplicates.bam' "
            "M=/dev/null >'/tmp/xxx/picard.errs' 2>&1",
            e.log[-2])
        self.assertEqual("$ mv '/tmp/xxx/picard-duplicates.bam' "
                         "'/tmp/xxx/result.bam'", e.log[-1])

    @patch('os.path.exists')
    def testMarkDuplicatesGATKBowtie2Threads(self, existsMock):
        """
        Using GATK to mark duplicates must result in the expected commands
        being run and output being produced. The number of threads is taken
        from the Bowtie2 instance.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp, threads=4)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq')
        bt.makeBAM()
        bt.indexBAM()
        bt.markDuplicatesGATK()
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nMarking duplicates with GATK.\n'))
        self.assertEqual(
            "$ gatk MarkDuplicatesSpark -I '/tmp/xxx/result.bam' -O "
            "'/tmp/xxx/gatk-duplicates.bam' --conf spark.executor.cores=4",
            e.log[-2])
        self.assertEqual("$ mv '/tmp/xxx/gatk-duplicates.bam' "
                         "'/tmp/xxx/result.bam'", e.log[-1])

    @patch('os.path.exists')
    def testMarkDuplicatesGATKThreadsInArg(self, existsMock):
        """
        Using GATK to mark duplicates must result in the expected commands
        being run and output being produced. The number of threads is given
        in the call to markDuplicatesGATK and must override the number given
        to the Bowtie2 instance.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp, threads=4)
        bt.buildIndex('file.fasta')
        bt.align(fastq1='file1.fastq')
        bt.makeBAM()
        bt.indexBAM()
        bt.markDuplicatesGATK(threads=32)
        log = fp.getvalue()
        self.assertTrue(log.endswith('\nMarking duplicates with GATK.\n'))
        self.assertEqual(
            "$ gatk MarkDuplicatesSpark -I '/tmp/xxx/result.bam' -O "
            "'/tmp/xxx/gatk-duplicates.bam' --conf spark.executor.cores=32",
            e.log[-2])
        self.assertEqual("$ mv '/tmp/xxx/gatk-duplicates.bam' "
                         "'/tmp/xxx/result.bam'", e.log[-1])

    def testClose(self):
        """
        Calling close() must result in the temporary directory being removed.
        """
        e = Executor(dryRun=True)
        fp = StringIO()
        bt = Bowtie2(executor=e, dryRun=True, verboseFp=fp)
        bt.close()
        self.assertEqual("$ rm -r '/tmp/xxx'", e.log[-1])
