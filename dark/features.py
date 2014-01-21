import numpy as np
import matplotlib.pyplot as plt

from dark.entrez import getSequence


class _Feature(object):
    """
    An offset-adjusted feature, with methods to return a textual description
    and a legend label.

    @param feature: A BioPython feature.
    @param offsetAdjuster: a function for adjusting feature X axis offsets for
        plotting.
    @param subfeature: A C{bool} to indicate if a feature is actually a
        subfeature.
    """

    def __init__(self, feature, offsetAdjuster, subfeature=False):
        self.feature = feature
        self.color = None  # Should be set with setColor
        self._offsetAdjuster = offsetAdjuster
        self.subfeature = subfeature

    def start(self):
        """
        Return the offset-adjusted start location of the feature.
        """
        return self._offsetAdjuster(int(self.feature.location.start))

    def end(self):
        """
        Return the offset-adjusted end location of the feature.
        """
        return self._offsetAdjuster(int(self.feature.location.end))

    def setColor(self, color):
        """
        An explicit method to set a feature's (plotting) color.

        @param color: A C{str} color.
        """
        self.color = color

    def description(self, excludedQualifiers=None, maxValueLength=None):
        """
        Provide a textual description of the feature and its qualifiers.

        @param excludedQualifiers: A C{set} of qualifier names that should be
            left out of the description. If C{None}, all feature qualifiers
            will be included.
        @param maxValueLength: An C{int} length, beyond which qualifier values
            will be truncated.
        @return: A C{str} description of the feature. The start and end offsets
            in the description are not offset-adjusted because the offset-
            adjusted values do not correspond to anything meaningful.
        """
        result = []
        excludedQualifiers = excludedQualifiers or set()
        for qualifier in sorted(self.feature.qualifiers.keys()):
            if qualifier not in excludedQualifiers:
                value = ', '.join(self.feature.qualifiers[qualifier])
                if qualifier == 'site_type' and value == 'other':
                    continue
                if maxValueLength and len(value) > maxValueLength:
                    value = value[:maxValueLength - 3] + '...'
                result.append('%s: %s' % (qualifier, value))
        return '%d-%d %s%s: %s' % (int(self.feature.location.start),
                                   int(self.feature.location.end),
                                   self.feature.type,
                                   ' (subfeature)' if self.subfeature else '',
                                   ', '.join(result))

    def legendLabel(self):
        return self.description(set(('db_xref', 'region_name')),
                                maxValueLength=30)


class _FeatureList(list):
    """
    Provide access to a list of L{Feature} objects.

    @param title: A C{str} sequence title from a BLAST hit. Of the form
        'gi|63148399|gb|DQ011818.1| Description...'.
    @param database: The S{str} name of the Entrez database to search.
    @param wantedTypes: A C{tuple} of feature types that are of interest.
        Feature whose types are not in this list will be ignored.
    @param offsetAdjuster: a function for adjusting feature X axis offsets for
        plotting.
    """

    def __init__(self, title, database, wantedTypes, offsetAdjuster):
        list.__init__(self)
        record = getSequence(title, db=database)
        if record is None:
            self.offline = True
        else:
            self.offline = False
            wantedTypes = set(wantedTypes)
            for feature in record.features:
                if feature.type in wantedTypes and feature.qualifiers:
                    self.append(_Feature(feature, offsetAdjuster))
                for subfeature in feature.sub_features:
                    if (subfeature.type in wantedTypes and
                            subfeature.qualifiers):
                        self.append(_Feature(subfeature, offsetAdjuster,
                                             subfeature=True))

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

    def add(self, fig, title, minX, maxX, offsetAdjuster):
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
        @return: If we seem to be offline, return C{None}. Otherwise, return
            a L{_FeatureList} instance.
        """

        fig.set_title('Target sequence features', fontsize=self.TITLE_FONTSIZE)
        fig.set_yticks([])

        features = _FeatureList(title, self.DATABASE, self.WANTED_TYPES,
                                offsetAdjuster)

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
            self._displayFeatures(fig, features, minX, maxX)
        else:
            self.tooManyFeaturesToPlot = True
            # fig.text(minX + (maxX - minX) / 3.0, 0,
            # 'Too many features to plot.', fontsize=self.FONTSIZE)
            fig.text(0.5, 0.5, 'Too many features to plot',
                     horizontalalignment='center', verticalalignment='center',
                     fontsize=self.FONTSIZE, transform=fig.transAxes)
            fig.axis([minX, maxX, -1, 1])

        return features

    def _displayFeatures(self, fig, features, minX, maxX):
        """
        Add the given C{features} to the figure in C{fig}.

        @param fig: A matplotlib figure.
        @param features: A C{_FeatureList} instance.
        @param minX: The smallest x coordinate.
        @param maxX: The largest x coordinate.
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

    def _displayFeatures(self, fig, features, minX, maxX):
        """
        Add the given C{features} to the figure in C{fig}.

        @param fig: A matplotlib figure.
        @param features: A C{_FeatureList} instance.
        @param minX: The smallest x coordinate.
        @param maxX: The largest x coordinate.
        """
        labels = []
        for index, feature in enumerate(features):
            fig.plot([feature.start(), feature.end()],
                     [index * -0.2, index * -0.2], color=feature.color,
                     linewidth=2)
            labels.append(feature.legendLabel())

        fig.axis([minX, maxX, (len(features) + 1) * -0.2, 0.2])

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
    WANTED_TYPES = ('CDS', 'mat_peptide', 'rRNA')

    def _displayFeatures(self, fig, features, minX, maxX):
        """
        Add the given C{features} to the figure in C{fig}.

        @param fig: A matplotlib figure.
        @param features: A C{_FeatureList} instance.
        @param minX: The smallest x coordinate.
        @param maxX: The largest x coordinate.
        """
        frame = None
        labels = []
        for feature in features:
            start = feature.start()
            end = feature.end()
            gene = feature.feature.qualifiers.get('gene', ['<no gene>'])[0]
            product = feature.feature.qualifiers.get(
                'product', ['<no product>'])[0]
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
                if product.lower().find('polyprotein') > -1:
                    y = frame + 0.2
                else:
                    y = frame
            fig.plot([start, end], [y, y], color=feature.color, linewidth=2)
            labels.append('%d-%d: %s (%s)' % (start, end, gene, product))

        fig.axis([minX, maxX, -1, 6])
        fig.set_yticks(np.arange(3))
        fig.set_ylabel('Frame', fontsize=17)
        if labels:
            # fig.legend(labels, bbox_to_anchor=(0.0, 1.1, 1.0, 0.102), loc=3,
            # ncol=3, mode='expand', borderaxespad=0.)
            fig.legend(labels, loc='upper left', ncol=3, shadow=True)
