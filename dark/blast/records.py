from dark.blast.hsp import printHSP


def printBlastRecord(record):
    """
    Print a BLAST record.

    @param record: A BioPython C{Bio.Blast.Record.Blast} instance.
    """
    for key in sorted(record.__dict__):
        if key not in ["alignments", "descriptions", "reference"]:
            print("%s: %r" % (key, record.__dict__[key]))
    print("alignments: (%d in total):" % len(record.alignments))
    for i, alignment in enumerate(record.alignments):
        print("  description %d:" % (i + 1))
        for attr in ["accession", "bits", "e", "num_alignments", "score"]:
            print("    %s: %s" % (attr, getattr(record.descriptions[i], attr)))
        print("  alignment %d:" % (i + 1))
        for attr in "accession", "hit_def", "hit_id", "length", "title":
            print("    %s: %s" % (attr, getattr(alignment, attr)))
        print("    HSPs (%d in total):" % len(alignment.hsps))
        for hspIndex, hsp in enumerate(alignment.hsps, start=1):
            print("      hsp %d:" % hspIndex)
            printHSP(hsp, "        ")
