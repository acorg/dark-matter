import numpy as np

try:
    from matplotlib import pyplot as plt
except ImportError:
    import platform
    if platform.python_implementation() == 'PyPy':
        # PyPy doesn't have a version of matplotlib. Make fake classes and
        # a Line2D function and that raise if used. This allows us to use
        # other 'dark' code that happens to import dark.mutations but not
        # use the functions that rely on matplotlib.
        class plt(object):
            def __getattr__(self, _):
                raise NotImplementedError(
                    'matplotlib is not supported under pypy')
    else:
        raise

from dark.entrez import getSequence


class Feature(object):
    """
    An offset-adjusted feature, with start and stop attributes and methods to
    return a textual description and a legend label.

    @param feature: A BioPython feature.
    @param subfeature: A C{bool} to indicate if a feature is actually a
        subfeature.
    """

    def __init__(self, feature, subfeature=False):
        self.feature = feature
        self.color = None  # Should be set with setColor
        self.subfeature = subfeature
        self.start = int(feature.location.start)
        self.end = int(feature.location.end)

    def setColor(self, color):
        """
        An explicit method to set a feature's (plotting) color.

        @param color: A C{str} color.
        """
        self.color = color

    def legendLabel(self):
        """
        Provide a textual description of the feature and its qualifiers to be
        used as a label in a plot legend.

        @return: A C{str} description of the feature.
        """
        excludedQualifiers = set((
            'codon_start', 'db_xref', 'protein_id', 'region_name',
            'ribosomal_slippage', 'rpt_type', 'translation', 'transl_except',
            'transl_table')
        )
        maxValueLength = 30
        result = []
        if self.feature.qualifiers:
            for qualifier in sorted(self.feature.qualifiers.keys()):
                if qualifier not in excludedQualifiers:
                    value = ', '.join(self.feature.qualifiers[qualifier])
                    if qualifier == 'site_type' and value == 'other':
                        continue
                    if len(value) > maxValueLength:
                        value = value[:maxValueLength - 3] + '...'
                    result.append('%s: %s' % (qualifier, value))
        return '%d-%d %s%s.%s' % (
            int(self.feature.location.start),
            int(self.feature.location.end),
            self.feature.type,
            ' (subfeature)' if self.subfeature else '',
            ' ' + ', '.join(result) if result else '')


class FeatureList(list):
    """
    Provide access to a list of L{Feature} objects.

    @param title: A C{str} sequence title from a BLAST hit. Of the form
        'gi|63148399|gb|DQ011818.1| Description...'.
    @param database: The S{str} name of the Entrez database to search.
    @param wantedTypes: A C{tuple} of feature types that are of interest.
        Feature whose types are not in this list will be ignored.
    @param sequenceFetcher: A function that takes a sequence title and a
        database name and returns a C{Bio.SeqIO} instance. If C{None}, use
        L{dark.entrez.getSequence}.
    """

    def __init__(self, title, database, wantedTypes, sequenceFetcher=None):
        list.__init__(self)
        self.offline = False
        sequenceFetcher = sequenceFetcher or getSequence
        try:
            record = sequenceFetcher(title, db=database)
        except ValueError:
            # Ignore. See https://github.com/acorg/dark-matter/issues/124
            return
        if record is None:
            self.offline = True
        else:
            wantedTypes = set(wantedTypes)
            for feature in record.features:
                if feature.type in wantedTypes:
                    self.append(Feature(feature))
                for subfeature in feature.sub_features:
                    if subfeature.type in wantedTypes:
                        self.append(Feature(subfeature, subfeature=True))

            # Assign colors to features.
            colormap = plt.cm.coolwarm
            colors = [colormap(i) for i in np.linspace(0.0, 0.99, len(self))]
            for feature, color in zip(self, colors):
                feature.setColor(color)


class _FeatureAdder(object):
    """
    Look up features for a title, and provide a method to add them to a figure
    as well as returning them.
    """

    TITLE_FONTSIZE = 16
    FONTSIZE = 20
    MAX_FEATURES_TO_DISPLAY = 50
    DATABASE = None  # Set in subclasses.
    WANTED_TYPES = None  # Set in subclasses.

    def __init__(self):
        self.tooManyFeaturesToPlot = False

    def add(self, fig, title, minX, maxX, offsetAdjuster=None,
            sequenceFetcher=None):
        """
        Find the features for a sequence title. If there aren't too many, add
        the features to C{fig}. Return information about the features, as
        described below.

        @param fig: A matplotlib figure.
        @param title: A C{str} sequence title from a BLAST hit. Of the form
            'gi|63148399|gb|DQ011818.1| Description...'.
        @param minX: The smallest x coordinate.
        @param maxX: The largest x coordinate.
         @param offsetAdjuster: a function for adjusting feature X axis offsets
            for plotting.
        @param sequenceFetcher: A function that takes a sequence title and a
            database name and returns a C{Bio.SeqIO} instance. If C{None}, use
            L{dark.entrez.getSequence}.
        @return: If we seem to be offline, return C{None}. Otherwise, return
            a L{FeatureList} instance.
        """

        offsetAdjuster = offsetAdjuster or (lambda x: x)

        fig.set_title('Target sequence features', fontsize=self.TITLE_FONTSIZE)
        fig.set_yticks([])

        features = FeatureList(title, self.DATABASE, self.WANTED_TYPES,
                               sequenceFetcher=sequenceFetcher)

        if features.offline:
            fig.text(minX + (maxX - minX) / 3.0, 0,
                     'You (or Genbank) appear to be offline.',
                     fontsize=self.FONTSIZE)
            fig.axis([minX, maxX, -1, 1])
            return None

        # If no interesting features were found, display a message saying
        # so in the figure.  Otherwise, if we don't have too many features
        # to plot, add the feature info to the figure.
        nFeatures = len(features)
        if nFeatures == 0:
            # fig.text(minX + (maxX - minX) / 3.0, 0, 'No features found',
            #          fontsize=self.FONTSIZE)
            fig.text(0.5, 0.5, 'No features found',
                     horizontalalignment='center', verticalalignment='center',
                     transform=fig.transAxes, fontsize=self.FONTSIZE)
            fig.axis([minX, maxX, -1, 1])
        elif nFeatures <= self.MAX_FEATURES_TO_DISPLAY:
            # Call the method in our subclass to do the figure display.
            self._displayFeatures(fig, features, minX, maxX, offsetAdjuster)
        else:
            self.tooManyFeaturesToPlot = True
            # fig.text(minX + (maxX - minX) / 3.0, 0,
            # 'Too many features to plot.', fontsize=self.FONTSIZE)
            fig.text(0.5, 0.5, 'Too many features to plot',
                     horizontalalignment='center', verticalalignment='center',
                     fontsize=self.FONTSIZE, transform=fig.transAxes)
            fig.axis([minX, maxX, -1, 1])

        return features

    def _displayFeatures(self, fig, features, minX, maxX, offsetAdjuster):
        """
        Add the given C{features} to the figure in C{fig}.

        @param fig: A matplotlib figure.
        @param features: A C{FeatureList} instance.
        @param minX: The smallest x coordinate.
        @param maxX: The largest x coordinate.
        @param offsetAdjuster: a function for adjusting feature X axis offsets
            for plotting.
        """
        raise NotImplementedError('_displayFeatures must be implemented in '
                                  'a subclass.')


class ProteinFeatureAdder(_FeatureAdder):
    """
    Subclass L{_FeatureAdder} with a method to add protein features to a
    figure.
    """
    DATABASE = 'protein'
    WANTED_TYPES = ('CDS', 'mat_peptide', 'rRNA', 'Site', 'Region')

    def _displayFeatures(self, fig, features, minX, maxX, offsetAdjuster):
        """
        Add the given C{features} to the figure in C{fig}.

        @param fig: A matplotlib figure.
        @param features: A C{FeatureList} instance.
        @param minX: The smallest x coordinate.
        @param maxX: The largest x coordinate.
        @param offsetAdjuster: a function for adjusting feature X axis offsets
            for plotting.
        """
        labels = []
        for index, feature in enumerate(features):
            fig.plot([offsetAdjuster(feature.start),
                      offsetAdjuster(feature.end)],
                     [index * -0.2, index * -0.2], color=feature.color,
                     linewidth=2)
            labels.append(feature.legendLabel())

        # Note that minX and maxX do not need to be adjusted by the offset
        # adjuster. They are the already-adjusted min/max values as
        # computed in computePlotInfo in blast.py
        fig.axis([minX, maxX, (len(features) + 1) * -0.2, 0.2])

        if labels:
            # Put a legend above the figure.
            box = fig.get_position()
            fig.set_position([box.x0, box.y0,
                              box.width, box.height * 0.2])
            fig.legend(labels, loc='lower center', bbox_to_anchor=(0.5, 1.4),
                       fancybox=True, shadow=True, ncol=2)


class NucleotideFeatureAdder(_FeatureAdder):
    """
    Subclass L{_FeatureAdder} with a method to add nucleotide features to a
    figure.
    """

    DATABASE = 'nucleotide'
    WANTED_TYPES = ('CDS', 'LTR', 'mat_peptide', 'misc_feature',
                    'misc_structure', 'repeat_region', 'rRNA')

    def _displayFeatures(self, fig, features, minX, maxX, offsetAdjuster):
        """
        Add the given C{features} to the figure in C{fig}.

        @param fig: A matplotlib figure.
        @param features: A C{FeatureList} instance.
        @param minX: The smallest x coordinate.
        @param maxX: The largest x coordinate.
        @param offsetAdjuster: a function for adjusting feature X axis offsets
            for plotting.
        """
        frame = None
        labels = []
        for feature in features:
            start = offsetAdjuster(feature.start)
            end = offsetAdjuster(feature.end)
            if feature.subfeature:
                subfeatureFrame = start % 3
                if subfeatureFrame == frame:
                    # Move overlapping subfeatures down a little to make them
                    # visible.
                    y = subfeatureFrame - 0.2
                else:
                    y = subfeatureFrame
            else:
                frame = start % 3
                # If we have a polyprotein, shift it up slightly so we can see
                # its components below it.
                product = feature.feature.qualifiers.get('product', [''])[0]
                if product.lower().find('polyprotein') > -1:
                    y = frame + 0.2
                else:
                    y = frame
            fig.plot([start, end], [y, y], color=feature.color, linewidth=2)
            labels.append(feature.legendLabel())

        # Note that minX and maxX do not need to be adjusted by the offset
        # adjuster. They are the already-adjusted min/max values as
        # computed in computePlotInfo in blast.py
        fig.axis([minX, maxX, -0.5, 2.5])
        fig.set_yticks(np.arange(3))
        fig.set_ylabel('Frame')

        if labels:
            # Put a legend above the figure.
            box = fig.get_position()
            fig.set_position([box.x0, box.y0,
                              box.width, box.height * 0.3])
            fig.legend(labels, loc='lower center', bbox_to_anchor=(0.5, 2.5),
                       fancybox=True, shadow=True, ncol=2)
