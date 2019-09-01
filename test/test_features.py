import numpy as np
from unittest import TestCase, skipUnless

try:
    import matplotlib
except ImportError:
    import platform
    if platform.python_implementation() == 'PyPy':
        havePlt = False
        # PyPy doesn't have a version of matplotlib. Make a fake class that
        # raises if used.

        class plt(object):
            def __getattr__(self, _):
                raise NotImplementedError(
                    'matplotlib is not supported under pypy')
    else:
        raise
else:
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    havePlt = True

try:
    from unittest.mock import call, MagicMock, ANY
except ImportError:
    from mock import call, MagicMock, ANY

from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from random import uniform

from dark.features import (Feature, FeatureList, _FeatureAdder,
                           ProteinFeatureAdder, NucleotideFeatureAdder)


def _randomLabel():
    return str(uniform(1, 1e9))


class Test_Feature(TestCase):
    """
    Tests of the C{Feature} class.
    """

    def testNotSubfeatureByDefault(self):
        """
        A feature must not be a subfeature by default.
        """
        location = FeatureLocation(100, 200)
        seqFeature = SeqFeature(location=location)
        feature = Feature(seqFeature)
        self.assertFalse(feature.subfeature)

    def testSubfeature(self):
        """
        A feature must be a subfeature if it is passed C{subfeature} = C{True}
        on initialization.
        """
        location = FeatureLocation(100, 200)
        seqFeature = SeqFeature(location=location)
        feature = Feature(seqFeature, subfeature=True)
        self.assertTrue(feature.subfeature)

    def testStart(self):
        """
        The start attribute must hold the feature's start.
        """
        location = FeatureLocation(100, 200)
        seqFeature = SeqFeature(location=location)
        feature = Feature(seqFeature)
        self.assertEqual(100, feature.start)

    def testEnd(self):
        """
        The end attribute must hold the feature start.
        """
        location = FeatureLocation(100, 200)
        seqFeature = SeqFeature(location=location)
        feature = Feature(seqFeature)
        self.assertEqual(200, feature.end)

    def testSetColor(self):
        """
        The setColor method must set the 'color' attribute on the feature.
        """
        location = FeatureLocation(100, 200)
        seqFeature = SeqFeature(location=location)
        feature = Feature(seqFeature)
        feature.setColor('red')
        self.assertEqual('red', feature.color)

    def testLegendLabel(self):
        """
        The legendLabel method must return a description.
        """
        location = FeatureLocation(100, 200)
        qualifiers = {
            'note': ['Capsid protein'],
        }
        seqFeature = SeqFeature(location=location, type='site',
                                qualifiers=qualifiers)
        feature = Feature(seqFeature)
        self.assertEqual('100-200 site. note: Capsid protein',
                         feature.legendLabel())

    def testLegendLabelQualifiersSorted(self):
        """
        The legendLabel method must return a description with sorted
        qualifiers.
        """
        location = FeatureLocation(100, 200)
        qualifiers = {
            'note': ['Capsid protein'],
            'product': ['CP1'],
        }
        seqFeature = SeqFeature(location=location, type='site',
                                qualifiers=qualifiers)
        feature = Feature(seqFeature)
        self.assertEqual('100-200 site. note: Capsid protein, product: CP1',
                         feature.legendLabel())

    def testLegendLabelNoQualifiers(self):
        """
        The legendLabel method must return a correct description for a feature
        that has no qualifiers.
        """
        location = FeatureLocation(100, 200)
        seqFeature = SeqFeature(location=location, type='site')
        feature = Feature(seqFeature)
        self.assertEqual('100-200 site.', feature.legendLabel())

    def testLegendLabelTruncatesValues(self):
        """
        The legendLabel method must return a description with qualifier
        values truncated to length 30 (including the trailing ...).
        """
        location = FeatureLocation(100, 200)
        qualifiers = {
            'note': ['x' * 40],
        }
        seqFeature = SeqFeature(location=location, type='site',
                                qualifiers=qualifiers)
        feature = Feature(seqFeature)
        xs = 'x' * 27 + '...'
        self.assertEqual('100-200 site. note: %s' % xs,
                         feature.legendLabel())


# We can't test anything that uses a FeatureList under pypy because
# dark.features makes use of matplotlib.
@skipUnless(havePlt, 'matplotlib not supported under pypy')
class Test_FeatureList(TestCase):
    """
    Tests of the C{FeatureList} class.
    """

    def testOffline(self):
        """
        If the sequence fetcher returns C{None} we must be marked as being
        offline.
        """
        def fetcher(title, db='database'):
            return None

        featureList = FeatureList('title', 'database', set(),
                                  sequenceFetcher=fetcher)
        self.assertEqual(True, featureList.offline)

    def testOfflineLength(self):
        """
        If the sequence fetcher indicates we're offline, the feature list
        must have length zero.
        """
        def fetcher(title, db='database'):
            return None

        featureList = FeatureList('title', 'database', set(),
                                  sequenceFetcher=fetcher)
        self.assertEqual(0, len(featureList))

    def testFetcherValueError(self):
        """
        If the sequence fetcher throws a C{ValueError}, the feature list
        must have length zero and should not be marked as being offline.
        """
        def fetcher(title, db):
            raise ValueError()

        featureList = FeatureList('title', 'database', set(),
                                  sequenceFetcher=fetcher)
        self.assertEqual(0, len(featureList))
        self.assertEqual(False, featureList.offline)

    def testNotOffline(self):
        """
        If the sequence fetcher does not return C{None} we must be marked as
        being online.
        """
        def fetcher(title, db='database'):
            return SeqRecord(None)

        featureList = FeatureList('title', 'database', set(),
                                  sequenceFetcher=fetcher)
        self.assertEqual(False, featureList.offline)

    def testNoFeaturesLength(self):
        """
        If the sequence fetcher returns a record with no features, the
        L{FeatureList} instance must have length zero.
        """
        def fetcher(title, db='database'):
            return SeqRecord(None)

        featureList = FeatureList('title', 'database', set(),
                                  sequenceFetcher=fetcher)
        self.assertEqual(0, len(featureList))

    def testNoQualifiersLength(self):
        """
        If the sequence fetcher returns a record with two features but
        neither of them has any qualifiers, the L{FeatureList} instance
        must still include both features.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='site', location=location)
            return SeqRecord(None, features=[feature, feature])

        featureList = FeatureList('title', 'database', set(['site']),
                                  sequenceFetcher=fetcher)
        self.assertEqual(2, len(featureList))

    def testNotSubfeature(self):
        """
        If the sequence fetcher returns a record with a feature that is
        not a subfeature, it must not be marked as a subfeature.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='site', qualifiers={'a': ['b']},
                                 location=location)
            return SeqRecord(None, features=[feature])

        featureList = FeatureList('title', 'database', set(['site']),
                                  sequenceFetcher=fetcher)
        self.assertFalse(featureList[0].subfeature)

    def testWantedTypeLength(self):
        """
        If the sequence fetcher returns a record with two features but
        only one of them has a wanted type ('site'), the L{FeatureList}
        instance must have length one.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature1 = SeqFeature(type='region', location=location)
            feature2 = SeqFeature(type='site', qualifiers={'a': ['b']},
                                  location=location)
            return SeqRecord(None, features=[feature1, feature2])

        featureList = FeatureList('title', 'database', set(['site']),
                                  sequenceFetcher=fetcher)
        self.assertEqual(1, len(featureList))

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testColors(self):
        """
        If the sequence fetcher returns a record with 3 features, each
        must be assigned a correct color.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='site', qualifiers={'a': ['b']},
                                 location=location)
            return SeqRecord(None, features=[feature] * 3)

        featureList = FeatureList('title', 'database', set(['site']),
                                  sequenceFetcher=fetcher)
        colormap = plt.cm.coolwarm
        colors = [colormap(i) for i in np.linspace(0.0, 0.99, 3)]

        for i in range(3):
            self.assertEqual(colors[i], featureList[i].color)


class Test_FeatureAdder(TestCase):
    """
    Tests of the C{_FeatureAdder} class.
    """

    def testNotTooManyFeaturesByDefault(self):
        """
        A new instance of L{_FeatureAdder} must be marked on creation as not
        having too many features to plot.
        """
        featureAdder = _FeatureAdder()
        self.assertFalse(featureAdder.tooManyFeaturesToPlot)

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testTitle(self):
        """
        A L{_FeatureAdder} must set the title of its figure.
        """
        def fetcher(title, db='database'):
            return None

        featureAdder = _FeatureAdder()
        fig = plt.subplot(111, label=_randomLabel())
        fig.set_title = MagicMock()
        featureAdder.add(fig, 'title', 0, 100, sequenceFetcher=fetcher)
        fig.set_title.assert_called_with('Target sequence features',
                                         fontsize=16)

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testYTicks(self):
        """
        A L{_FeatureAdder} must set the title of its figure.
        """
        def fetcher(title, db='database'):
            return None

        featureAdder = _FeatureAdder()
        fig = plt.subplot(111, label=_randomLabel())
        fig.set_yticks = MagicMock()
        featureAdder.add(fig, 'title', 0, 100, sequenceFetcher=fetcher)
        fig.set_yticks.assert_called_with([])

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testOffline(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns C{None},
        the C{add} method of the L{_FeatureAdder} must return C{None}. The
        C{text} and C{axis} methods on the fig must both be called correctly.
        """
        def fetcher(title, db='database'):
            return None

        featureAdder = _FeatureAdder()
        fig = plt.subplot(111, label=_randomLabel())
        fig.text = MagicMock()
        fig.axis = MagicMock()
        featureAdder.add(fig, 'title', 0, 300, sequenceFetcher=fetcher)
        fig.text.assert_called_with(100, 0,
                                    'You (or Genbank) appear to be offline.',
                                    fontsize=20)
        fig.axis.assert_called_with([0, 300, -1, 1])

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testNoFeatures(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns no features
        the C{text} and C{axis} methods on the fig must both be called
        correctly and the C{add} method must return C{[]}.
        """
        def fetcher(title, db='database'):
            return SeqRecord(None)

        featureAdder = _FeatureAdder()
        featureAdder.WANTED_TYPES = ('site',)
        fig = plt.subplot(111, label=_randomLabel())
        fig.text = MagicMock()
        fig.axis = MagicMock()
        result = featureAdder.add(fig, 'title', 0, 300,
                                  sequenceFetcher=fetcher)
        fig.text.assert_called_with(0.5, 0.5,
                                    'No features found',
                                    horizontalalignment='center',
                                    verticalalignment='center',
                                    transform=ANY, fontsize=20)
        fig.axis.assert_called_with([0, 300, -1, 1])
        self.assertEqual([], result)

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testOneFeatureRaisesNotImplementedError(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns a feature,
        the C{add} method must raise C{NotImplementedError} because the
        C{_displayFeatures} method is not implemented.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='site', qualifiers={'a': ['b']},
                                 location=location)
            return SeqRecord(None, features=[feature])

        featureAdder = _FeatureAdder()
        featureAdder.WANTED_TYPES = ('site',)
        fig = plt.subplot(111, label=_randomLabel())
        self.assertRaises(NotImplementedError, featureAdder.add, fig, 'title',
                          0, 300, sequenceFetcher=fetcher)

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testTooManyFeatures(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns too many
        features, the C{text} and C{axis} methods on the figure must be called
        correctly and the C{add} call must return the sequences.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='site', qualifiers={'a': ['b']},
                                 location=location)
            return SeqRecord(None, features=[feature] * 100)

        featureAdder = _FeatureAdder()
        featureAdder.WANTED_TYPES = ('site',)
        fig = plt.subplot(111, label=_randomLabel())
        fig.text = MagicMock()
        fig.axis = MagicMock()
        result = featureAdder.add(fig, 'title', 0, 300,
                                  sequenceFetcher=fetcher)
        fig.text.assert_called_with(0.5, 0.5,
                                    'Too many features to plot',
                                    horizontalalignment='center',
                                    verticalalignment='center',
                                    transform=ANY, fontsize=20)
        fig.axis.assert_called_with([0, 300, -1, 1])
        self.assertTrue(isinstance(result, FeatureList))
        self.assertEqual(100, len(result))


class TestProteinFeatureAdder(TestCase):
    """
    Tests of the C{ProteinFeatureAdder} class.
    """

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testUnwantedFeature(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns a feature
        whose type is not wanted, the figure's plot method must not be called
        and the C{add} method must return an empty feature list.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='unwanted', qualifiers={'a': ['b']},
                                 location=location)
            return SeqRecord(None, features=[feature])

        featureAdder = ProteinFeatureAdder()
        fig = plt.subplot(111, label=_randomLabel())
        fig.plot = MagicMock()
        result = featureAdder.add(fig, 'title', 0, 300,
                                  sequenceFetcher=fetcher)
        self.assertEqual([], fig.plot.call_args_list)
        self.assertEqual([], result)

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testOneFeature(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns a feature,
        the C{text} and C{axis} methods on the figure must be called correctly
        and the C{add} call must return the sequences.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='Site', qualifiers={'a': ['b']},
                                 location=location)
            return SeqRecord(None, features=[feature])

        featureAdder = ProteinFeatureAdder()
        fig = plt.subplot(111, label=_randomLabel())
        fig.plot = MagicMock()
        fig.axis = MagicMock()
        fig.legend = MagicMock()
        result = featureAdder.add(fig, 'title', 0, 300,
                                  sequenceFetcher=fetcher)
        fig.plot.assert_called_with(
            [100, 200], [-0.0, -0.0],
            color=(0.2298057, 0.298717966, 0.75368315299999999, 1.0),
            linewidth=2)
        fig.axis.assert_called_with([0, 300, -0.4, 0.2])
        fig.legend.assert_called_with(
            ['100-200 Site. a: b'], loc='lower center', shadow=True,
            bbox_to_anchor=(0.5, 1.4), ncol=2, fancybox=True)
        self.assertTrue(isinstance(result, FeatureList))
        self.assertEqual(1, len(result))


class TestNucleotideFeatureAdder(TestCase):
    """
    Tests of the C{NucleotideFeatureAdder} class.
    """

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testUnwantedFeature(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns a feature
        whose type is not wanted, the figure's plot method must not be called
        and the C{add} method must return an empty feature list.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='unwanted', qualifiers={'a': ['b']},
                                 location=location)
            return SeqRecord(None, features=[feature])

        featureAdder = NucleotideFeatureAdder()
        fig = plt.subplot(111, label=_randomLabel())
        fig.plot = MagicMock()
        result = featureAdder.add(fig, 'title', 0, 300,
                                  sequenceFetcher=fetcher)
        self.assertEqual([], fig.plot.call_args_list)
        self.assertEqual([], result)

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testOneFeature(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns a feature,
        the C{text} and C{axis} methods on the figure must be called correctly
        and the C{add} call must return the sequences.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='CDS', qualifiers={'a': ['b']},
                                 location=location)
            return SeqRecord(None, features=[feature])

        featureAdder = NucleotideFeatureAdder()
        fig = plt.subplot(111, label=_randomLabel())
        fig.plot = MagicMock()
        fig.axis = MagicMock()
        fig.legend = MagicMock()
        result = featureAdder.add(fig, 'title', 0, 300,
                                  sequenceFetcher=fetcher)
        fig.plot.assert_called_with(
            [100, 200], [1, 1],
            color=(0.2298057, 0.298717966, 0.75368315299999999, 1.0),
            linewidth=2)
        fig.axis.assert_called_with([0, 300, -0.5, 2.5])
        fig.legend.assert_called_with(
            ['100-200 CDS. a: b'], loc='lower center',
            shadow=True, ncol=2, fancybox=True,
            bbox_to_anchor=(0.5, 2.5))
        self.assertTrue(isinstance(result, FeatureList))
        self.assertEqual(1, len(result))

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testOneFeatureAdjusted(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns a feature,
        the C{text} and C{axis} methods on the figure must be called correctly
        and the C{add} call must return the sequences.

        Note that offsets in the legend of nucleotide plots are adjusted. They
        shouldn't be as the adjusted offsets make no sense to the reader.
        """
        def fetcher(title, db='database'):
            location = FeatureLocation(100, 200)
            feature = SeqFeature(type='CDS', qualifiers={'a': ['b']},
                                 location=location)
            return SeqRecord(None, features=[feature])

        featureAdder = NucleotideFeatureAdder()
        fig = plt.subplot(111, label=_randomLabel())
        fig.plot = MagicMock()
        fig.axis = MagicMock()
        fig.legend = MagicMock()

        def adjuster(x):
            return 3 * x

        result = featureAdder.add(fig, 'title', 0, 300, adjuster,
                                  sequenceFetcher=fetcher)
        fig.plot.assert_called_with(
            [300, 600], [0, 0],
            color=(0.2298057, 0.298717966, 0.75368315299999999, 1.0),
            linewidth=2)
        fig.axis.assert_called_with([0, 300, -0.5, 2.5])
        fig.legend.assert_called_with(
            ['100-200 CDS. a: b'], loc='lower center',
            shadow=True, ncol=2, fancybox=True,
            bbox_to_anchor=(0.5, 2.5))
        self.assertTrue(isinstance(result, FeatureList))
        self.assertEqual(1, len(result))

    @skipUnless(havePlt, 'matplotlib not supported under pypy')
    def testPolyproteinsAreMovedUp(self):
        """
        If the sequence fetcher used by a L{_FeatureAdder} returns a feature,
        that's a polyprotein, the feature must be plotted a little (0.2) above
        its normal location.
        """
        def fetcher(title, db='database'):
            feature1 = SeqFeature(type='CDS',
                                  qualifiers={'product': ['a polyprotein']},
                                  location=FeatureLocation(100, 200))
            feature2 = SeqFeature(type='CDS',
                                  qualifiers={'a': ['b']},
                                  location=FeatureLocation(130, 150))
            return SeqRecord(None, features=[feature1, feature2])

        featureAdder = NucleotideFeatureAdder()
        fig = plt.subplot(111, label=_randomLabel())
        fig.plot = MagicMock()
        fig.axis = MagicMock()
        fig.legend = MagicMock()
        result = featureAdder.add(fig, 'title', 0, 300,
                                  sequenceFetcher=fetcher)
        self.assertEqual(
            fig.plot.call_args_list,
            [
                call([100, 200], [1.2, 1.2], color=ANY, linewidth=2),
                call([130, 150], [1.0, 1.0], color=ANY, linewidth=2),
            ])
        fig.axis.assert_called_with([0, 300, -0.5, 2.5])
        fig.legend.assert_called_with(
            ['100-200 CDS. product: a polyprotein',
             '130-150 CDS. a: b'],
            loc='lower center', shadow=True, ncol=2, fancybox=True,
            bbox_to_anchor=(0.5, 2.5))
        self.assertTrue(isinstance(result, FeatureList))
        self.assertEqual(2, len(result))
