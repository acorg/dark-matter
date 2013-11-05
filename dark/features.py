import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

START_CODONS = set(['ATG'])
STOP_CODONS = set(['TAA', 'TAG', 'TGA'])


def findCodons(seq, codons):
    """
    Find all instances of the codons in 'codons' in the given sequence.

    seq: A Bio.Seq.Seq instance.
    codons: A set of codon strings.

    Return: a generator yielding matching codon offsets.
    """
    seqLen = len(seq)
    start = 0
    while start < seqLen:
        triplet = str(seq[start:start + 3])
        if triplet in codons:
            yield start
        start = start + 3


def getFeatures(fig, record, minX, maxX):
    """
    fig is a matplotlib figure.
    record is a Bio.Seq with features, or None (if offline).
    minX: the smallest x coordinate.
    maxX: the largest x coordinate.
    """
    fig.set_title('Target sequence features', fontsize=20)
    # print record.features

    toPlot = []
    totalSubfeatures = 0
    if record:
        for feature in record.features:
            if feature.type in ('CDS', 'rRNA'):
                toPlot.append(feature)
                totalSubfeatures += len(feature.sub_features)

    if len(toPlot) > 20:
        fig.text(minX + (maxX - minX) / 3.0, 0,
                 ('Too many features to plot.'),
                 fontsize=16)
        fig.axis([minX, maxX, -1, 1])
        fig.set_yticks([])
        return toPlot

    elif record is None or not toPlot:
        fig.text(minX + (maxX - minX) / 3.0, 0,
                 ('No features found.' if record
                  else 'You (or Genbank) appear to be offline.'),
                 fontsize=16)
        fig.axis([minX, maxX, -1, 1])
        fig.set_yticks([])
        return []

    else:
        return toPlot


def addFeatures(fig, record, minX, maxX):

    toPlot = getFeatures(fig, record, minX, maxX)
    result = []
    plotFeatures = True

    if len(toPlot) > 20:
        plotFeatures = False
        for feature in toPlot:
            featureType = str(feature.type)
            start = int(feature.location.start)
            end = int(feature.location.end)
            gene = feature.qualifiers.get('gene', ['<no gene>'])[0]
            product = feature.qualifiers.get('product', ['<no product>'])[0]
            result.append('%d-%d: %s (%s)' % (start, end, gene, product))
            for subfeature in feature.sub_features:
                start = int(subfeature.location.start)
                end = int(subfeature.location.end)
                subfeatureFrame = start % 3
                result.append('%d-%d: %s subfeature' % (start, end, gene))
        return plotFeatures, result

    else:
        colormap = plt.cm.coolwarm
        colors = [colormap(i) for i in
                  np.linspace(0.0, 0.99, len(toPlot) + totalSubfeatures)]
        labels = []

        index = -1
        for feature in toPlot:
            index += 1
            start = int(feature.location.start)
            end = int(feature.location.end)
            result.append({
                'color': colors[index],
                'end': end,
                'start': start,
            })
            frame = start % 3
            fig.plot([start, end], [frame, frame], color=colors[index],
                     linewidth=2)
            gene = feature.qualifiers.get('gene', ['<no gene>'])[0]
            product = feature.qualifiers.get('product', ['<no product>'])[0]
            labels.append('%d-%d: %s (%s)' % (start, end, gene, product))
            for subfeature in feature.sub_features:
                index += 1
                start = int(subfeature.location.start)
                end = int(subfeature.location.end)
                result.append({
                    'color': colors[index],
                    'end': end,
                    'start': start,
                })
                subfeatureFrame = start % 3
                if subfeatureFrame == frame:
                    # Move overlapping subfeatures down a little to make them
                    # visible.
                    subfeatureFrame -= 0.2
                fig.plot([start, end], [subfeatureFrame, subfeatureFrame],
                         color=colors[index])
                labels.append('%d-%d: %s subfeature' % (start, end, gene))

        fig.axis([minX, maxX, -1, 6])
        fig.set_yticks(np.arange(3))
        fig.set_ylabel('Frame', fontsize=17)
        if labels:
            # fig.legend(labels, bbox_to_anchor=(0.0, 1.1, 1.0, 0.102), loc=3,
            # ncol=3, mode='expand', borderaxespad=0.)
            fig.legend(labels, loc='upper left', ncol=3, shadow=True)

    return plotFeatures, result


def addORFs(fig, seq, minX, maxX, featureEndpoints):
    """
    fig is a matplotlib figure.
    seq is a Bio.Seq.Seq.
    minX: the smallest x coordinate.
    maxX: the largest x coordinate.
    featureEndpoints: an array of features as returned by addFeatures (may be
        empty).
    """
    for frame in range(3):
        target = seq[frame:]
        for (codons, codonType, color) in (
                (START_CODONS, 'start', 'green'),
                (STOP_CODONS, 'stop', 'red')):
            offsets = list(findCodons(target, codons))
            if offsets:
                fig.plot(offsets, np.tile(frame, len(offsets)), marker='.',
                         markersize=4, color=color, linestyle='None')

    # Add the feature endpoints.
    for fe in featureEndpoints:
        line = Line2D([fe['start'], fe['start']], [-1, 3], color=fe['color'],
                      linewidth=1)
        fig.add_line(line)
        line = Line2D([fe['end'], fe['end']], [-1, 3], color='#cccccc')
        fig.add_line(line)

    fig.axis([minX, maxX, -1, 3])
    fig.set_yticks(np.arange(3))
    fig.set_ylabel('Frame', fontsize=17)
    fig.set_title('Target sequence start (%s) and stop (%s) codons' % (
        ', '.join(sorted(START_CODONS)), ', '.join(sorted(STOP_CODONS))),
        fontsize=20)


def addReversedORFs(fig, seq, minX, maxX):
    """
    fig is a matplotlib figure.
    seq is a Bio.Seq.Seq (the reverse complement of the sequence we're
        plotting against).
    minX: the smallest x coordinate.
    maxX: the largest x coordinate.
    """
    for frame in range(3):
        target = seq[frame:]
        for (codons, codonType, color) in (
                (START_CODONS, 'start', 'green'),
                (STOP_CODONS, 'stop', 'red')):
            offsets = map(lambda offset: maxX - offset,
                          findCodons(target, codons))
            if offsets:
                fig.plot(offsets, np.tile(frame, len(offsets)), marker='.',
                         markersize=4, color=color, linestyle='None')

    fig.axis([minX, maxX, -1, 3])
    fig.set_yticks(np.arange(3))
    fig.set_ylabel('Frame', fontsize=17)
    fig.set_title('Reversed target sequence start (%s) & stop (%s) codons' % (
        ', '.join(sorted(START_CODONS)), ', '.join(sorted(STOP_CODONS))),
        fontsize=20)
