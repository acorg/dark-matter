import numpy as np

START_CODONS = set(["ATG"])
STOP_CODONS = set(["TAA", "TAG", "TGA"])


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
        triplet = str(seq[start : start + 3])
        if triplet in codons:
            yield start
        start = start + 3


def addORFs(fig, seq, minX, maxX, offsetAdjuster):
    """
    fig is a matplotlib figure.
    seq is a Bio.Seq.Seq.
    minX: the smallest x coordinate.
    maxX: the largest x coordinate.
    featureEndpoints: an array of features as returned by addFeatures (may be
        empty).
    offsetAdjuster: a function to adjust feature X axis offsets for plotting.
    """
    for frame in range(3):
        target = seq[frame:]
        for codons, codonType, color in (
            (START_CODONS, "start", "green"),
            (STOP_CODONS, "stop", "red"),
        ):
            offsets = list(map(offsetAdjuster, findCodons(target, codons)))
            if offsets:
                fig.plot(
                    offsets,
                    np.tile(frame, len(offsets)),
                    marker=".",
                    markersize=4,
                    color=color,
                    linestyle="None",
                )

    fig.axis([minX, maxX, -1, 3])
    fig.set_yticks(np.arange(3))
    fig.set_ylabel("Frame", fontsize=17)
    fig.set_title(
        "Subject start (%s) and stop (%s) codons"
        % (", ".join(sorted(START_CODONS)), ", ".join(sorted(STOP_CODONS))),
        fontsize=20,
    )


def addReversedORFs(fig, seq, minX, maxX, offsetAdjuster):
    """
    fig is a matplotlib figure.
    seq is a Bio.Seq.Seq (the reverse complement of the sequence we're
        plotting against).
    minX: the smallest x coordinate.
    maxX: the largest x coordinate.
    offsetAdjuster: a function to adjust feature X axis offsets for plotting.
    """
    for frame in range(3):
        target = seq[frame:]
        for codons, codonType, color in (
            (START_CODONS, "start", "green"),
            (STOP_CODONS, "stop", "red"),
        ):
            offsets = [
                maxX - offsetAdjuster(offset) for offset in findCodons(target, codons)
            ]
            if offsets:
                fig.plot(
                    offsets,
                    np.tile(frame, len(offsets)),
                    marker=".",
                    markersize=4,
                    color=color,
                    linestyle="None",
                )

    fig.axis([minX, maxX, -1, 3])
    fig.set_yticks(np.arange(3))
    fig.set_ylabel("Frame", fontsize=17)
    fig.set_title(
        "Reversed subject start (%s) and stop (%s) codons"
        % (", ".join(sorted(START_CODONS)), ", ".join(sorted(STOP_CODONS))),
        fontsize=20,
    )
