from unittest import TestCase
from six import StringIO, assertRaisesRegex
from six.moves import builtins
from contextlib import contextmanager

from dark.proteins import (
    splitNames, _NO_PATHOGEN_NAME, getPathogenProteinCounts, ProteinGrouper,
    PathogenSampleFiles)

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from .mocking import File


class TestSplitNames(TestCase):
    """
    Tests for the splitNames function.
    """

    def testNoBrackets(self):
        """
        If a string with no trailing pathogen name in square brackets is given,
        splitNames must return its argument followed by a string indicating
        no pathogen name could be found.
        """
        self.assertEqual(('xxx', _NO_PATHOGEN_NAME), splitNames('xxx'))

    def testTwoSetsOfBrackets(self):
        """
        If a string with two trailing substrings in square brackets is given,
        splitNames must extract the substring from the second set of square
        brackets and use that as the pathogen name.
        """
        self.assertEqual(('xxx [other]', 'pathogen'),
                         splitNames('xxx [other] [pathogen]'))

    def testWhitespaceStripping(self):
        """
        If a string with names that have whitespace is passed, splitNames
        must strip the whitespace in its result.
        """
        self.assertEqual(('xxx', 'pathogen'),
                         splitNames(' xxx [ pathogen ]'))

    def testNestedBrackets(self):
        """
        If a string with two nested trailing substrings in square brackets is
        given, splitNames must return its argument followed by a string
        indicating no pathogen name could be found.
        """
        self.assertEqual(('xxx [nested [pathogen name]]', _NO_PATHOGEN_NAME),
                         splitNames('xxx [nested [pathogen name]]'))

    def testNormalCase(self):
        """
        If a string with a protein and pathogen name is passed, splitNames
        must return the expected result.
        """
        self.assertEqual(('protein name', 'pathogen name'),
                         splitNames('protein name [pathogen name]'))


class TestGetPathogenProteinCounts(TestCase):
    """
    Tests for the getPathogenProteinCounts function.
    """
    def testNone(self):
        """
        getPathogenProteinCounts must return an empty result if passed None as
        the protein FASTA file.
        """
        self.assertEqual({}, getPathogenProteinCounts(None))

    def testExpected(self):
        """
        getPathogenProteinCounts must return the expected result.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename.fasta', filename)
                    self.count += 1
                    return File(['>protein 1 [pathogen 1]\n',
                                 'ACTG\n',
                                 '>protein 2 [pathogen 1]\n',
                                 'AA\n',
                                 '>no pathogen name here\n',
                                 'AA\n',
                                 '>protein 3 [pathogen 2]\n',
                                 'AA\n'])
                else:
                    self.test.fail('We are only supposed to be called once!')

        sideEffect = SideEffect(self)
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            self.assertEqual(
                {
                    'pathogen 1': 2,
                    'pathogen 2': 1,
                },
                getPathogenProteinCounts(['filename.fasta']))
            self.assertEqual(1, sideEffect.count)

    def testExpectedWithTwoFiles(self):
        """
        getPathogenProteinCounts must return the expected result when details
        are read from two FASTA files.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('filename1.fasta', filename)
                    self.count += 1
                    return File(['>protein 1 [pathogen 1]\n',
                                 'ACTG\n',
                                 '>protein 3 [pathogen 2]\n',
                                 'AA\n'])
                elif self.count == 1:
                    self.test.assertEqual('filename2.fasta', filename)
                    self.count += 1
                    return File(['>protein 2 [pathogen 1]\n',
                                 'AA\n'])
                else:
                    self.test.fail('We are only supposed to be called twice!')

        sideEffect = SideEffect(self)
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            self.assertEqual(
                {
                    'pathogen 1': 2,
                    'pathogen 2': 1,
                },
                getPathogenProteinCounts(
                    ['filename1.fasta', 'filename2.fasta']))
            self.assertEqual(2, sideEffect.count)


class TestProteinGrouper(TestCase):
    """
    Tests for the dark.proteins.ProteinGrouper class.
    """

    def testUnknownFormat(self):
        """
        Passing an unknown format argument must result in a ValueError
        being raised.
        """
        error = "^format_ must be either 'fasta' or 'fastq'\.$"
        assertRaisesRegex(self, ValueError, error, ProteinGrouper,
                          format_='unknown')

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
        pathogen names attributes must both be empty.
        """
        pg = ProteinGrouper()
        self.assertEqual({}, pg.pathogenNames)
        self.assertEqual({}, pg.sampleNames)

    def testUnknownPathogenType(self):
        """
        If the toHTML method of a protein grouper is given an unknown pathogen
        type it must raise a ValueError.
        """
        pg = ProteinGrouper()
        error = ("^Unrecognized pathogenType argument: 'x'\. Value must be "
                 "either 'bacterial' or 'viral'\.$")
        assertRaisesRegex(self, ValueError, error, pg.toHTML, pathogenType='x')

    def testDuplicatePathogenProteinSample(self):
        """
        If a protein grouper is given duplicate information for a
        pathogen/protein/sample combination it must raise a ValueError.
        """
        fp = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327|X|I44.6 ubiquitin [Lausannevirus]\n')
        pg = ProteinGrouper()
        pg.addFile('sample', fp)
        fp.seek(0)
        error = ("^Protein 'gi\|327\|X\|I44.6 ubiquitin' already seen for "
                 "pathogen 'Lausannevirus' sample 'sample'\.$")
        assertRaisesRegex(self, ValueError, error, pg.addFile, 'sample', fp)

    def testOneLineInOneFile(self):
        """
        If a protein grouper is given one file with one line, its pathogenNames
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
                        'proteins': {
                            'gi|327|X|I44.6 ubiquitin': {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.77,
                                'readsFilename': 'out/0.fasta',
                                'hspCount': 6,
                                'index': 0,
                                'medianScore': 46.6,
                                'outDir': 'out',
                                'proteinLength': 74,
                                'proteinName': 'gi|327|X|I44.6 ubiquitin',
                                'proteinURL': (
                                    'http://www.ncbi.nlm.nih.gov/nuccore/I44'),
                                'readCount': 5,
                            },
                        },
                        'uniqueReadCount': None,
                    },
                }
            },
            pg.pathogenNames)

    def testOneLineInOneFileFASTQ(self):
        """
        If a protein grouper is given one file with one line, its pathogenNames
        dict must be as expected, including for a FASTQ file.
        """
        fp = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327|X|I44.6 ubiquitin [Lausannevirus]\n')
        pg = ProteinGrouper(format_='fastq')
        pg.addFile('sample-filename', fp)
        self.assertEqual(
            {
                'Lausannevirus': {
                    'sample-filename': {
                        'proteins': {
                            'gi|327|X|I44.6 ubiquitin': {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.77,
                                'readsFilename': 'out/0.fastq',
                                'hspCount': 6,
                                'index': 0,
                                'medianScore': 46.6,
                                'outDir': 'out',
                                'proteinLength': 74,
                                'proteinName': 'gi|327|X|I44.6 ubiquitin',
                                'proteinURL': (
                                    'http://www.ncbi.nlm.nih.gov/nuccore/I44'),
                                'readCount': 5,
                            },
                        },
                        'uniqueReadCount': None,
                    },
                }
            },
            pg.pathogenNames)

    def testOneLineInOneFileTitle(self):
        """
        If a protein grouper is given one file with one line, its _title method
        must return the expected string.
        """
        fp = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327|X|I44.6 ubiquitin [Lausannevirus]\n')
        pg = ProteinGrouper()
        pg.addFile('sample-filename', fp)
        self.assertEqual(
            'Overall, proteins from 1 pathogen were found in 1 sample.',
            pg._title())

    def testTwoLinesInOneFileTitle(self):
        """
        If a protein grouper is given one file with two protein lines, each
        from a different pathogen, its _title method must return the expected
        string.
        """
        fp = StringIO(
            '0.77 46.6 48.1 5 6 74 gi|327|X|I44.6 ubiquitin [Lausannevirus]\n'
            '0.77 46.6 48.1 5 6 74 gi|327|X|I44.6 ubiquitin [X Virus]\n'
            )
        pg = ProteinGrouper()
        pg.addFile('sample-filename', fp)
        self.assertEqual(
            'Overall, proteins from 2 pathogens were found in 1 sample.',
            pg._title())

    def testTwoLinesInOneFileSamePathogen(self):
        """
        If a protein grouper is given one file with two lines from the same
        pathogen, its pathogenNames dict must be as expected.
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
                        'proteins': {
                            'gi|327410| protein 77': {
                                'bestScore': 44.2,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.63,
                                'readsFilename': 'out/0.fasta',
                                'hspCount': 9,
                                'index': 0,
                                'medianScore': 41.3,
                                'outDir': 'out',
                                'proteinLength': 12,
                                'proteinName': 'gi|327410| protein 77',
                                'proteinURL': None,
                                'readCount': 9,
                            },
                            'gi|327409| ubiquitin': {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'out/1.png',
                                'coverage': 0.77,
                                'readsFilename': 'out/1.fasta',
                                'hspCount': 6,
                                'index': 1,
                                'medianScore': 46.6,
                                'outDir': 'out',
                                'proteinLength': 74,
                                'proteinName': 'gi|327409| ubiquitin',
                                'proteinURL': None,
                                'readCount': 5,
                            },
                        },
                        'uniqueReadCount': None,
                    },
                },
            },
            pg.pathogenNames)

    def testTwoLinesInOneFileDifferentPathogens(self):
        """
        If a protein grouper is given one file with two lines from different
        pathogens, its pathogenNames dict must be as expected.
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
                        'proteins': {
                            'gi|327410| protein 77': {
                                'bestScore': 44.2,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.63,
                                'readsFilename': 'out/0.fasta',
                                'hspCount': 9,
                                'index': 0,
                                'medianScore': 41.3,
                                'outDir': 'out',
                                'proteinLength': 12,
                                'proteinName': 'gi|327410| protein 77',
                                'proteinURL': None,
                                'readCount': 9,
                            },
                        },
                        'uniqueReadCount': None,
                    },
                },
                'Hepatitis B virus': {
                    'sample-filename': {
                        'proteins': {
                            'gi|327409| ubiquitin': {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'out/1.png',
                                'coverage': 0.77,
                                'readsFilename': 'out/1.fasta',
                                'hspCount': 6,
                                'index': 1,
                                'medianScore': 46.6,
                                'outDir': 'out',
                                'proteinLength': 74,
                                'proteinName': 'gi|327409| ubiquitin',
                                'proteinURL': None,
                                'readCount': 5,
                            },
                        },
                        'uniqueReadCount': None,
                    },
                },
            },
            pg.pathogenNames)

    def testOneLineInEachOfTwoFilesSamePathogen(self):
        """
        If a protein grouper is given two files, each with one line from the
        same pathogen, its pathogenNames dict must be as expected.
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
                        'proteins': {
                            'gi|327410| protein 77': {
                                'bestScore': 44.2,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.63,
                                'readsFilename': 'out/0.fasta',
                                'hspCount': 9,
                                'index': 0,
                                'medianScore': 41.3,
                                'outDir': 'out',
                                'proteinLength': 12,
                                'proteinName': 'gi|327410| protein 77',
                                'proteinURL': None,
                                'readCount': 9,
                            },
                        },
                        'uniqueReadCount': None,
                    },
                    'sample-filename-2': {
                        'proteins': {
                            'gi|327409| ubiquitin': {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'out/0.png',
                                'coverage': 0.77,
                                'readsFilename': 'out/0.fasta',
                                'hspCount': 6,
                                'index': 0,
                                'medianScore': 46.6,
                                'outDir': 'out',
                                'proteinLength': 74,
                                'proteinName': 'gi|327409| ubiquitin',
                                'proteinURL': None,
                                'readCount': 5,
                            },
                        },
                        'uniqueReadCount': None,
                    },
                },
            },
            pg.pathogenNames)

    def testOneLineInEachOfTwoFilesSamePathogenTitle(self):
        """
        If a protein grouper is given two files, each with one line from the
        same pathogen, its _title method must return the expected string.
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
            'Overall, proteins from 1 pathogen were found in 2 samples.',
            pg._title())

    def testOneLineInEachOfTwoFilesDifferentPathogens(self):
        """
        If a protein grouper is given two files in two different directories,
        each with one line from the different pathogens, its pathogenNames dict
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
                        'proteins': {
                            'gi|327410| protein 77': {
                                'bestScore': 44.2,
                                'bluePlotFilename': 'dir-1/out/0.png',
                                'coverage': 0.63,
                                'readsFilename': 'dir-1/out/0.fasta',
                                'hspCount': 9,
                                'index': 0,
                                'medianScore': 41.3,
                                'outDir': 'dir-1/out',
                                'proteinLength': 12,
                                'proteinName': 'gi|327410| protein 77',
                                'proteinURL': None,
                                'readCount': 9,
                            },
                        },
                        'uniqueReadCount': None,
                    },
                },
                'Hepatitis B virus': {
                    'dir-2/sample-filename-2': {
                        'proteins': {
                            'gi|327409| ubiquitin': {
                                'bestScore': 48.1,
                                'bluePlotFilename': 'dir-2/out/0.png',
                                'coverage': 0.77,
                                'readsFilename': 'dir-2/out/0.fasta',
                                'hspCount': 6,
                                'index': 0,
                                'medianScore': 46.6,
                                'outDir': 'dir-2/out',
                                'proteinLength': 74,
                                'proteinName': 'gi|327409| ubiquitin',
                                'proteinURL': None,
                                'readCount': 5,
                            },
                        },
                        'uniqueReadCount': None,
                    },
                },
            },
            pg.pathogenNames)

    def testOneLineInEachOfTwoFilesDifferentPathogensTitle(self):
        """
        If a protein grouper is given two files, each with one line from
        different pathogens, its _title method must return the expected string.
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
        self.assertEqual(
            'Overall, proteins from 2 pathogens were found in 2 samples.',
            pg._title())

    def testNoFilesToStr(self):
        """
        If no files have been given to a protein grouper, its text string
        format must as expected.
        """
        pg = ProteinGrouper()
        self.assertEqual(
            'Overall, proteins from 0 pathogens were found in 0 samples.\n',
            pg.toStr())

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
            'Overall, proteins from 1 pathogen were found in 1 sample.\n'
            '\n'
            'HBV (in 1 sample)\n'
            '  sample-filename (1 protein, 5 reads)\n'
            '    0.77\t46.60\t48.10\t   5\t   6\t  0\tgi|32|X|I4 protein X\n',
            pg.toStr())

    def testMaxProteinFraction(self):
        """
        The maxProteinFraction method must return the correct values.
        """
        class SideEffect(object):
            def __init__(self, test):
                self.test = test
                self.count = 0

            def sideEffect(self, filename, **kwargs):
                if self.count == 0:
                    self.test.assertEqual('proteins.fasta', filename)
                    self.count += 1
                    return File(['>protein 1 [pathogen 1]\n',
                                 'ACTG\n',
                                 '>protein 2 [pathogen 1]\n',
                                 'AA\n',
                                 '>protein 3 [pathogen 1]\n',
                                 'AA\n',
                                 '>protein 4 [pathogen 1]\n',
                                 'AA\n',
                                 '>no pathogen name here\n',
                                 'AA\n',
                                 '>protein 5 [pathogen 2]\n',
                                 'AA\n'])
                else:
                    self.test.fail('We are only supposed to be called once!')

        sideEffect = SideEffect(self)
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            pg = ProteinGrouper(proteinFastaFilenames=['proteins.fasta'])
            self.assertEqual(1, sideEffect.count)

            fp = StringIO(
                '0.77 46.6 48.1 5 6 74 gi|32|X|I4 protein 1 [pathogen 1]\n'
                '0.77 46.6 48.1 5 6 74 gi|32|X|I4 protein 5 [pathogen 2]\n'
            )
            pg.addFile('sample-1', fp)

            fp = StringIO(
                '0.77 46.6 48.1 5 6 74 gi|32|X|I4 protein 2 [pathogen 1]\n'
                '0.77 46.6 48.1 5 6 74 gi|32|X|I4 protein 3 [pathogen 1]\n'
            )
            pg.addFile('sample-1', fp)

            fp = StringIO(
                '0.77 46.6 48.1 5 6 74 gi|32|X|I4 protein 1 [pathogen 1]\n'
                '0.77 46.6 48.1 5 6 74 gi|32|X|I4 protein 2 [pathogen 1]\n'
            )
            pg.addFile('sample-2', fp)

            self.assertEqual(0.75, pg.maxProteinFraction('pathogen 1'))
            self.assertEqual(1.0, pg.maxProteinFraction('pathogen 2'))


class TestPathogenSampleFiles(TestCase):
    """
    Tests for the PathogenSampleFiles class.
    """

    def testUnknownFormat(self):
        """
        Passing an unknown format argument must result in a ValueError
        being raised.
        """
        pg = ProteinGrouper()
        error = "^format_ must be either 'fasta' or 'fastq'\.$"
        assertRaisesRegex(self, ValueError, error, PathogenSampleFiles,
                          pg, format_='unknown')

    def testOpenNotCalledOnRepeatedCall(self):
        """
        If a repeated call to pathogenSampleFiles.add is made with the same
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
                    self.test.assertEqual('out/pathogen-0-sample-0.fasta',
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
        pathogenSampleFiles = PathogenSampleFiles(pg)

        sideEffect = Open(self, manager()).sideEffect
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = sideEffect
            filename = pathogenSampleFiles.add('Lausannevirus', 'filename-1')
            self.assertEqual('out/pathogen-0-sample-0.fasta', filename)
            self.assertEqual('>id1\nACTG\n', fastaIO.getvalue())

            # Repeated call. The side effect open will fail if open is
            # called at this point.
            filename = pathogenSampleFiles.add('Lausannevirus', 'filename-1')
            self.assertEqual('out/pathogen-0-sample-0.fasta', filename)

    def testIdenticalReadsRemoved(self):
        """
        If two proteins in the same pathogen are matched by the same read, the
        de-duplicated FASTA for the pathogen must have only one copy of the
        duplicated read.
        """
        class Open(object):
            def __init__(self, test, manager):
                self.test = test
                self.manager = manager
                self.expectedFilenames = {'out/0.fasta', 'out/1.fasta',
                                          'out/pathogen-0-sample-0.fasta'}

            def sideEffect(self, filename, *args, **kwargs):
                try:
                    self.expectedFilenames.remove(filename)
                except KeyError:
                    self.test.fail(
                        'Open called with unexpected filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))
                else:
                    if filename == 'out/0.fasta':
                        return File(['>id1\n', 'ACTG\n'])
                    elif filename == 'out/1.fasta':
                        return File(['>id1\n', 'ACTG\n', '>id2\n', 'CAGT\n'])
                    else:
                        return self.manager

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
        pathogenSampleFiles = PathogenSampleFiles(pg)

        opener = Open(self, manager())
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = opener.sideEffect
            filename = pathogenSampleFiles.add('Lausannevirus', 'filename-1')

        self.assertEqual('out/pathogen-0-sample-0.fasta', filename)
        self.assertEqual('>id1\nACTG\n>id2\nCAGT\n', fastaIO.getvalue())
        # Make sure all expected filenames were seen by the mocked open.
        self.assertEqual(set(), opener.expectedFilenames)

    def testProteinsSavedCorrectly(self):
        """
        Information about proteins must be saved correctly in the
        ProteinGrouper for a given pathogen/sample combination.
        """
        class Open(object):
            def __init__(self, test, manager):
                self.test = test
                self.manager = manager
                self.expectedFilenames = {'out/0.fasta', 'out/1.fasta',
                                          'out/pathogen-0-sample-0.fasta'}

            def sideEffect(self, filename, *args, **kwargs):
                if filename in self.expectedFilenames:
                    if filename == 'out/0.fasta':
                        return File(['>id1\n', 'ACTG\n'])
                    elif filename == 'out/1.fasta':
                        return File(['>id2\n', 'AC\n', '>id3\n', 'CAGTTTT\n'])
                    else:
                        return self.manager
                else:
                    self.test.fail(
                        'Open called with unexpected filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        fp = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
            '0.77 46.6 48.1 5 6 74 gi|327409| ubiquitin [Lausannevirus]\n'
        )
        fastaIO = StringIO()

        @contextmanager
        def manager():
            yield fastaIO

        opener = Open(self, manager())
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = opener.sideEffect
            pg = ProteinGrouper()
            pg.addFile('filename-1', fp)
            pathogenSampleFiles = PathogenSampleFiles(pg)
            pathogenSampleFiles.add('Lausannevirus', 'filename-1')

        self.assertEqual(
            {
                'proteins': {
                    'gi|327409| ubiquitin': {
                        'bestScore': 48.1,
                        'bluePlotFilename': 'out/1.png',
                        'coverage': 0.77,
                        'hspCount': 6,
                        'index': 1,
                        'medianScore': 46.6,
                        'outDir': 'out',
                        'proteinLength': 74,
                        'proteinName': 'gi|327409| ubiquitin',
                        'proteinURL': None,
                        'readCount': 5,
                        'readsFilename': 'out/1.fasta',
                    },
                    'gi|327410| protein 77': {
                        'bestScore': 44.2,
                        'bluePlotFilename': 'out/0.png',
                        'coverage': 0.63,
                        'hspCount': 9,
                        'index': 0,
                        'medianScore': 41.3,
                        'outDir': 'out',
                        'proteinLength': 12,
                        'proteinName': 'gi|327410| protein 77',
                        'proteinURL': None,
                        'readCount': 9,
                        'readsFilename': 'out/0.fasta',
                    }
                },
                'uniqueReadCount': 3,
            },
            pg.pathogenNames['Lausannevirus']['filename-1'])

    def testReadLengthsAdded(self):
        """
        If saveReadLengths is True for a ProteinGrouper, read lengths must be
        saved for each protein.
        """
        class Open(object):
            def __init__(self, test, manager):
                self.test = test
                self.manager = manager
                self.expectedFilenames = {'out/0.fasta', 'out/1.fasta',
                                          'out/pathogen-0-sample-0.fasta'}

            def sideEffect(self, filename, *args, **kwargs):
                if filename in self.expectedFilenames:
                    if filename == 'out/0.fasta':
                        return File(['>id1\n', 'ACTG\n'])
                    elif filename == 'out/1.fasta':
                        return File(['>id2\n', 'AC\n', '>id3\n', 'CAGTTTT\n'])
                    else:
                        return self.manager
                else:
                    self.test.fail(
                        'Open called with unexpected filename: %r, Args: %r, '
                        'Keyword args: %r.' % (filename, args, kwargs))

        fp = StringIO(
            '0.63 41.3 44.2 9 9 12 gi|327410| protein 77 [Lausannevirus]\n'
            '0.77 46.6 48.1 5 6 74 gi|327409| ubiquitin [Lausannevirus]\n'
        )
        fastaIO = StringIO()

        @contextmanager
        def manager():
            yield fastaIO

        opener = Open(self, manager())
        with patch.object(builtins, 'open') as mockMethod:
            mockMethod.side_effect = opener.sideEffect
            pg = ProteinGrouper(saveReadLengths=True)
            pg.addFile('filename-1', fp)
            pathogenSampleFiles = PathogenSampleFiles(pg)
            pathogenSampleFiles.add('Lausannevirus', 'filename-1')

        # Read lengths must be saved correctly.
        proteins = pg.pathogenNames['Lausannevirus']['filename-1']['proteins']
        self.assertEqual((4,),
                         proteins['gi|327410| protein 77']['readLengths'])
        self.assertEqual((2, 7),
                         proteins['gi|327409| ubiquitin']['readLengths'])
