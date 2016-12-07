from unittest import TestCase
from six import StringIO
from six.moves import builtins
from contextlib import contextmanager

from dark.proteins import ProteinGrouper, VirusSampleFASTA

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from .mocking import File


class TestProteinGrouper(TestCase):
    """
    Tests for the dark.proteins.ProteinGrouper class
    """

    def testNoAssetDir(self):
        """
        If no asset directorey is given to a protein grouper, its _assetDir
        attribute be the default ('out').
        """
        pg = ProteinGrouper()
        self.assertEqual('out', pg._assetDir)

    def testAssetDir(self):
        """
        If an asset directorey is given to a protein grouper, its _assetDir
        attribute be set to hold that value.
        """
        pg = ProteinGrouper(assetDir='xxx')
        self.assertEqual('xxx', pg._assetDir)

    def testNoRegex(self):
        """
        If no regex is given to a protein grouper, its _sampleNameRegex
        attribute be None.
        """
        pg = ProteinGrouper()
        self.assertEqual(None, pg._sampleNameRegex)

    def testNoFiles(self):
        """
        If no files have been given to a protein grouper, its sample names and
        virus titles attributes must both be empty.
        """
        pg = ProteinGrouper()
        self.assertEqual({}, pg.virusTitles)
        self.assertEqual({}, pg.sampleNames)

    def testOneLineInOneFile(self):
        """
        If a protein grouper is given one file with one line, its virusTitles
        dict must be as expected.
        """
        fp = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327|X|I44.6 ubiquitin [Lausannevirus]\n')
        pg = ProteinGrouper()
        pg.addFile('sample-filename', fp)
        self.assertEqual(
            {
                'Lausannevirus': {
                    'sample-filename': {
                        'proteins': [
                            {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.77,
                                'fastaFilename': 'out/0.fasta',
                                'hspCount': 6,
                                'index': 0,
                                'medianScore': 46.6,
                                'outDir': 'out',
                                'proteinLength': 74,
                                'proteinTitle': 'gi|327|X|I44.6 ubiquitin',
                                'proteinURL': (
                                    'http://www.ncbi.nlm.nih.gov/nuccore/I44'),
                                'readCount': 5,
                            },
                        ],
                        'uniqueReadCount': None,
                    },
                }
            },
            pg.virusTitles)

    def testOneLineInOneFileTitle(self):
        """
        If a protein grouper is given one file with one line, its _title method
        must return the expected string.
        """
        fp = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327|X|I44.6 ubiquitin [Lausannevirus]\n')
        pg = ProteinGrouper()
        pg.addFile('sample-filename', fp)
        self.assertEqual('1 virus found in 1 sample', pg._title())

    def testTwoLinesInOneFileTitle(self):
        """
        If a protein grouper is given one file with two protein lines, each
        from a different virus, its _title method must return the expected
        string.
        """
        fp = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327|X|I44.6 ubiquitin [Lausannevirus]\n'
            '0.77 46.6 48.1 5 6 74 gi|327|X|I44.6 ubiquitin [X Virus]\n'
            )
        pg = ProteinGrouper()
        pg.addFile('sample-filename', fp)
        self.assertEqual('2 viruses found in 1 sample', pg._title())

    def testTwoLinesInOneFileSameVirus(self):
        """
        If a protein grouper is given one file with two lines from the same
        virus, its virusTitles dict must be as expected.
        """
        fp = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
            '0.77 46.6 48.1 5 6 74 gi|327409| ubiquitin [Lausannevirus]\n'
        )
        pg = ProteinGrouper()
        pg.addFile('sample-filename', fp)
        self.assertEqual(
            {
                'Lausannevirus': {
                    'sample-filename': {
                        'proteins': [
                            {
                                'bestScore': 44.2,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.63,
                                'fastaFilename': 'out/0.fasta',
                                'hspCount': 9,
                                'index': 0,
                                'medianScore': 41.3,
                                'outDir': 'out',
                                'proteinLength': 12,
                                'proteinTitle': 'gi|327410| protein 77',
                                'proteinURL': None,
                                'readCount': 9,
                            },
                            {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'out/1.png',
                                'coverage': 0.77,
                                'fastaFilename': 'out/1.fasta',
                                'hspCount': 6,
                                'index': 1,
                                'medianScore': 46.6,
                                'outDir': 'out',
                                'proteinLength': 74,
                                'proteinTitle': 'gi|327409| ubiquitin',
                                'proteinURL': None,
                                'readCount': 5,
                            },
                        ],
                        'uniqueReadCount': None,
                    },
                },
            },
            pg.virusTitles)

    def testTwoLinesInOneFileDifferentViruses(self):
        """
        If a protein grouper is given one file with two lines from different
        viruss, its virusTitles dict must be as expected.
        """
        fp = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
            '0.77 46.6 48.1 5 6 74 gi|327409| ubiquitin [Hepatitis B virus]\n'
        )
        pg = ProteinGrouper()
        pg.addFile('sample-filename', fp)
        self.assertEqual(
            {
                'Lausannevirus': {
                    'sample-filename': {
                        'proteins': [
                            {
                                'bestScore': 44.2,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.63,
                                'fastaFilename': 'out/0.fasta',
                                'hspCount': 9,
                                'index': 0,
                                'medianScore': 41.3,
                                'outDir': 'out',
                                'proteinLength': 12,
                                'proteinTitle': 'gi|327410| protein 77',
                                'proteinURL': None,
                                'readCount': 9,
                            },
                        ],
                        'uniqueReadCount': None,
                    },
                },
                'Hepatitis B virus': {
                    'sample-filename': {
                        'proteins': [
                            {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'out/1.png',
                                'coverage': 0.77,
                                'fastaFilename': 'out/1.fasta',
                                'hspCount': 6,
                                'index': 1,
                                'medianScore': 46.6,
                                'outDir': 'out',
                                'proteinLength': 74,
                                'proteinTitle': 'gi|327409| ubiquitin',
                                'proteinURL': None,
                                'readCount': 5,
                            },
                        ],
                        'uniqueReadCount': None,
                    },
                },
            },
            pg.virusTitles)

    def testOneLineInEachOfTwoFilesSameVirus(self):
        """
        If a protein grouper is given two files, each with one line from the
        same virus, its virusTitles dict must be as expected.
        """
        fp1 = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
        )
        fp2 = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327409| ubiquitin [Lausannevirus]\n'
        )
        pg = ProteinGrouper()
        pg.addFile('sample-filename-1', fp1)
        pg.addFile('sample-filename-2', fp2)
        self.assertEqual(
            {
                'Lausannevirus': {
                    'sample-filename-1': {
                        'proteins': [
                            {
                                'bestScore': 44.2,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.63,
                                'fastaFilename': 'out/0.fasta',
                                'hspCount': 9,
                                'index': 0,
                                'medianScore': 41.3,
                                'outDir': 'out',
                                'proteinLength': 12,
                                'proteinTitle': 'gi|327410| protein 77',
                                'proteinURL': None,
                                'readCount': 9,
                            },
                        ],
                        'uniqueReadCount': None,
                    },
                    'sample-filename-2': {
                        'proteins': [
                            {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.77,
                                'fastaFilename': 'out/0.fasta',
                                'hspCount': 6,
                                'index': 0,
                                'medianScore': 46.6,
                                'outDir': 'out',
                                'proteinLength': 74,
                                'proteinTitle': 'gi|327409| ubiquitin',
                                'proteinURL': None,
                                'readCount': 5,
                            },
                        ],
                        'uniqueReadCount': None,
                    },
                },
            },
            pg.virusTitles)

    def testOneLineInEachOfTwoFilesSameVirusTitle(self):
        """
        If a protein grouper is given two files, each with one line from the
        same virus, its _title method must return the expected string.
        """
        fp1 = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
        )
        fp2 = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327409| ubiquitin [Lausannevirus]\n'
        )
        pg = ProteinGrouper()
        pg.addFile('sample-filename-1', fp1)
        pg.addFile('sample-filename-2', fp2)
        self.assertEqual('1 virus found in 2 samples', pg._title())

    def testOneLineInEachOfTwoFilesDifferentViruses(self):
        """
        If a protein grouper is given two files in two different directories,
        each with one line from the different viruses, its virusTitles dict
        must be as expected.
        """
        fp1 = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
        )
        fp2 = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327409| ubiquitin [Hepatitis B virus]\n'
        )
        pg = ProteinGrouper()
        pg.addFile('dir-1/sample-filename-1', fp1)
        pg.addFile('dir-2/sample-filename-2', fp2)
        self.assertEqual(
            {
                'Lausannevirus': {
                    'dir-1/sample-filename-1': {
                        'proteins': [
                            {
                                'bestScore': 44.2,
                                'bluePlotFilename': 'dir-1/out/0.png',
                                'coverage': 0.63,
                                'fastaFilename': 'dir-1/out/0.fasta',
                                'hspCount': 9,
                                'index': 0,
                                'medianScore': 41.3,
                                'outDir': 'dir-1/out',
                                'proteinLength': 12,
                                'proteinTitle': 'gi|327410| protein 77',
                                'proteinURL': None,
                                'readCount': 9,
                            },
                        ],
                        'uniqueReadCount': None,
                    },
                },
                'Hepatitis B virus': {
                    'dir-2/sample-filename-2': {
                        'proteins': [
                            {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'dir-2/out/0.png',
                                'coverage': 0.77,
                                'fastaFilename': 'dir-2/out/0.fasta',
                                'hspCount': 6,
                                'index': 0,
                                'medianScore': 46.6,
                                'outDir': 'dir-2/out',
                                'proteinLength': 74,
                                'proteinTitle': 'gi|327409| ubiquitin',
                                'proteinURL': None,
                                'readCount': 5,
                            },
                        ],
                        'uniqueReadCount': None,
                    },
                },
            },
            pg.virusTitles)

    def testOneLineInEachOfTwoFilesDifferentVirusesTitle(self):
        """
        If a protein grouper is given two files, each with one line from
        different viruses, its _title method must return the expected string.
        """
        fp1 = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
        )
        fp2 = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327409| ubiquitin [HBV]\n'
        )
        pg = ProteinGrouper()
        pg.addFile('sample-filename-1', fp1)
        pg.addFile('sample-filename-2', fp2)
        self.assertEqual('2 viruses found in 2 samples', pg._title())

    def testNoFilesToStr(self):
        """
        If no files have been given to a protein grouper, its text string
        format must as expected.
        """
        pg = ProteinGrouper()
        self.assertEqual('0 viruses found in 0 samples\n', pg.toStr())

    def testNoFilesToHTML(self):
        """
        If no files have been given to a protein grouper, its HTML string
        format must as expected.
        """
        pg = ProteinGrouper()
        self.assertEqual(
            '\n'.join([
                '<html>',
                '<head>',
                '<title>',
                '0 viruses found in 0 samples',
                '</title>',
                '</head>',
                '<body>',
                '<style>',
                '            body {',
                '                margin-left: 2%;',
                '                margin-right: 2%;',
                '            }',
                '            .sample {',
                '                margin-bottom: 2px;',
                '            }',
                '            .sample-name {',
                '                color: red;',
                '            }',
                '            .index {',
                '                font-size: small;',
                '            }',
                '            .protein-title {',
                '                font-family: "Courier New", Courier, '
                'monospace;',
                '            }',
                '            .stats {',
                '                font-family: "Courier New", Courier, '
                'monospace;',
                '                white-space: pre;',
                '            }',
                '            .protein-list {',
                '                margin-top: 2px;',
                '            }',
                '</style>',
                '</head>',
                '<body>',
                '<h1>0 viruses found in 0 samples</h1>',
                '<h2>Virus index</h2>',
                '</p>',
                '<h2>Sample index</h2>',
                '</p>',
                '<h1>Viruses by sample</h1>',
                '<p>',
                'In the bullet point protein lists below, there are eight '
                'fields:',
                '<ol>',
                '<li>Coverage fraction.</li>',
                '<li>Median bit score.</li>',
                '<li>Best bit score.</li>',
                '<li>Read count.</li>',
                '<li>HSP count (a read can match a protein more than once).'
                '</li>',
                '<li>Protein length (in AAs).</li>',
                '<li>Index (just ignore this).</li>',
                '<li>Protein name.</li>',
                '</ol>',
                '</p>',
                '<h1>Samples by virus</h1>',
                '<p>',
                'In the bullet point protein lists below, there are eight '
                'fields:',
                '<ol>',
                '<li>Coverage fraction.</li>',
                '<li>Median bit score.</li>',
                '<li>Best bit score.</li>',
                '<li>Read count.</li>',
                '<li>HSP count (a read can match a protein more than once).'
                '</li>',
                '<li>Protein length (in AAs).</li>',
                '<li>Index (just ignore this).</li>',
                '<li>Protein name.</li>',
                '</ol>',
                '</p>',
                '</body>',
                '</html>',
                ]),
            pg.toHTML())

    def testOneLineInOneFileToStr(self):
        """
        If a protein grouper is given one file with one line, its toStr method
        must produce the expected result.
        """
        fp = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|32|X|I4 protein X [HBV]\n')
        pg = ProteinGrouper()
        pg.addFile('sample-filename', fp)
        self.assertEqual(
            '1 virus found in 1 sample\n'
            '\n'
            'HBV (in 1 sample)\n'
            '  sample-filename (1 protein, 5 reads)\n'
            '    0.77\t46.60\t48.10\t   5\t   6\t  0\tgi|32|X|I4 protein X\n',
            pg.toStr())


class TestVirusSampleFASTA(TestCase):
    """
    Tests for the VirusSampleFASTA class.
    """

    def testOpenNotCalledOnRepeatedCall(self):
        """
        If a repeated call to virusSampleFASTA.add is made with the same
        arguments, no file should be read because the original result value is
        cached.
        """
        class Open(object):
            def __init__(self, test, manager):
                self.test = test
                self.manager = manager
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('out/0.fasta', filename)
                    self.count += 1
                    return File(['>id1\n', 'ACTG\n'])
                elif self.count == 1:
                    self.test.assertEqual('out/virus-0-sample-0.fasta',
                                          filename)
                    self.count += 1
                    return self.manager
                else:
                    self.test.fail(
                        'We are only supposed to be called twice. '
                        'Filename: %r, Args: %r, Keyword args: %r.' %
                        (filename, args, kwargs))

        fp = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
        )
        fastaIO = StringIO()

        @contextmanager
        def manager():
            yield fastaIO

        pg = ProteinGrouper()
        pg.addFile('filename-1', fp)
        virusSampleFASTA = VirusSampleFASTA(pg)

        sideEffect = Open(self, manager()).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            filename = virusSampleFASTA.add('Lausannevirus', 'filename-1')
            self.assertEqual('out/virus-0-sample-0.fasta', filename)
            self.assertEqual('>id1\nACTG\n', fastaIO.getvalue())

            # Repeated call. The side effect open will fail if open is
            # called at this point.
            filename = virusSampleFASTA.add('Lausannevirus', 'filename-1')
            self.assertEqual('out/virus-0-sample-0.fasta', filename)

    def testIdenticalReadsRemoved(self):
        """
        If two proteins in the same virus are matched by the same read, the
        de-duplicated FASTA for the virus must have only one copy of the
        duplicated read.
        """
        class Open(object):
            def __init__(self, test, manager):
                self.test = test
                self.manager = manager
                self.count = 0

            def sideEffect(self, filename, *args, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('out/0.fasta', filename)
                    self.count += 1
                    return File(['>id1\n', 'ACTG\n'])
                elif self.count == 1:
                    self.test.assertEqual('out/1.fasta', filename)
                    self.count += 1
                    return File(['>id1\n', 'ACTG\n', '>id2\n', 'CAGT\n'])
                elif self.count == 2:
                    self.test.assertEqual('out/virus-0-sample-0.fasta',
                                          filename)
                    self.count += 1
                    return self.manager
                else:
                    self.test.fail(
                        'We are only supposed to be called 3 times. '
                        'Filename: %r, Args: %r, Keyword args: %r.' %
                        (filename, args, kwargs))

        fp = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
            '0.77 46.6 48.1 5 6 74 gi|327409| ubiquitin [Lausannevirus]\n'
        )
        fastaIO = StringIO()

        @contextmanager
        def manager():
            yield fastaIO

        pg = ProteinGrouper()
        pg.addFile('filename-1', fp)
        virusSampleFASTA = VirusSampleFASTA(pg)

        sideEffect = Open(self, manager()).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            filename = virusSampleFASTA.add('Lausannevirus', 'filename-1')
            self.assertEqual('out/virus-0-sample-0.fasta', filename)
            self.assertEqual('>id1\nACTG\n>id2\nCAGT\n', fastaIO.getvalue())
