#!/usr/bin/env python

from typing import Dict

import sys
import argparse
from collections import OrderedDict, defaultdict
from operator import itemgetter

from dark.aligners import edlibAlign, mafft, needle
from dark.dna import compareDNAReads
from dark.filter import (
    addFASTAFilteringCommandLineOptions, parseFASTAFilteringCommandLineOptions,
    addFASTAEditingCommandLineOptions, parseFASTAEditingCommandLineOptions)
from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions

MAFFT_DEFAULT_ARGS = '--globalpair --maxiterate 1000 --preservecase'
MAFFT_ALGORITHMS_URL = (
    'https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html')
NEEDLE_DEFAULT_ARGS = 'auto'


def align(reads, args):
    """
    Align a pair of reads.

    @param reads: An iterable of two C{DNARead} instances.
    @param args: An argparse C{Namespace} instance with command-line options.
    @return: A C{list} of two aligned C{DNARead} instances.
    """
    if args.aligner == 'mafft':
        # Be careful in examining args.alignerOptions because we want the
        # user to be able to pass an empty string (so check against None
        # before deciding to use the default.)
        options = (MAFFT_DEFAULT_ARGS if args.alignerOptions is None
                   else args.alignerOptions)
        return mafft(reads, args.verbose, options=options,
                     threads=args.threads)
    elif args.aligner == 'needle':
        # Be careful in examining args.alignerOptions because we want the
        # user to be able to pass an empty string (so check against None
        # before deciding to use the default.)
        options = (NEEDLE_DEFAULT_ARGS if args.alignerOptions is None
                   else args.alignerOptions)
        return needle(reads, args.verbose, options=options)
    else:
        assert args.aligner == 'edlib'
        return edlibAlign(reads)


def thresholdToCssName(threshold):
    """
    Turn a floating point threshold into a string that can be used as a CSS
    class name.

    param threshold: The C{float} threshold.
    return: A C{str} CSS class name.
    """
    # '.' is illegal in a CSS class name.
    return 'threshold-' + str(threshold).replace('.', '_')


def thresholdForIdentity(identity, colors):
    """
    Get the best identity threshold for a specific identity value.

    @param identity: A C{float} nucleotide identity.
    @param colors: A C{list} of (threshold, color) tuples, where threshold is a
        C{float} and color is a C{str} to be used as a cell background. This
        is as returned by C{parseColors}.
    @return: The first C{float} threshold that the given identity is at least
        as big as.
    """
    for threshold, _ in colors:
        if identity >= threshold:
            return threshold
    raise ValueError('This should never happen! Last threshold is not 0.0?')


def parseColors(colors, defaultColor):
    """
    Parse command line color information.

    @param colors: A C{list} of space separated "value color" strings, such as
        ["0.9 red", "0.75 rgb(23, 190, 207)", "0.1 #CF3CF3"].
    @param defaultColor: The C{str} color to use for cells that do not reach
        the identity fraction threshold of any color in C{colors}.
    @return: A C{list} of (threshold, color) tuples, where threshold is a
        C{float} (from C{colors}) and color is a C{str} (from C{colors}). The
        list will be sorted by decreasing threshold values.
    """
    result = []
    if colors:
        for colorInfo in colors:
            fields = colorInfo.split(None, 1)
            if len(fields) == 2:
                threshold, color = fields
                try:
                    threshold = float(threshold)
                except ValueError:
                    print('--color arguments must be given as space-separated '
                          'pairs of "value color" where the value is a '
                          'numeric identity threshold. Your value %r is not '
                          'numeric.' % threshold, file=sys.stderr)
                    sys.exit(1)
                if 0.0 > threshold > 1.0:
                    print('--color arguments must be given as space-separated '
                          'pairs of "value color" where the value is a '
                          'numeric identity threshold from 0.0 to 1.0. Your '
                          'value %r is not in that range.' % threshold,
                          file=sys.stderr)
                    sys.exit(1)

                result.append((threshold, color))
            else:
                print('--color arguments must be given as space-separated '
                      'pairs of "value color". You have given %r, which does '
                      'not contain a space.' % colorInfo, file=sys.stderr)
                sys.exit(1)

    result.sort(key=itemgetter(0), reverse=True)

    if not result or result[-1][0] > 0.0:
        result.append((0.0, defaultColor))

    return result


def getReadLengths(reads, gapChars):
    """
    Get all read lengths, excluding gap characters.

    @param reads: A C{Reads} instance.
    @param gapChars: A C{str} of sequence characters considered to be gaps.
    @return: A C{dict} keyed by read id, with C{int} length values.
    """
    gapChars = set(gapChars)
    result = {}
    for read in reads:
        result[read.id] = len(read) - sum(
            character in gapChars for character in read.sequence)
    return result


def explanation(matchAmbiguous, concise, showLengths, showGaps, showNs):
    """
    Make an explanation of the output HTML table.

    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct. Otherwise, we are strict
        and insist that only non-ambiguous nucleotides can contribute to the
        matching nucleotide count.
    @param concise: If C{True}, do not show match detail abbreviations.
    @param showLengths: If C{True}, include the lengths of sequences.
    @param showGaps: If C{True}, include the number of gaps in sequences.
    @param showNs: If C{True}, include the number of N characters in sequences.
    @return: A C{str} of HTML.
    """
    result = ["""
<h1>Sequence versus sequence identity table</h1>

<p>

The table cells below show the nucleotide identity fraction for the sequences
(<span class="best">like this</span> for the best value in each row). The
identity fraction numerator is the sum of the number of identical
    """]

    if matchAmbiguous:
        result.append('nucleotides plus the number of ambiguously matching '
                      'nucleotides.')
    else:
        result.append('nucleotides.')

    result.append("""The denominator
is the length of the sequence <em>for the row</em>. Sequence gaps
are not included when calculating their lengths.

</p>
    """)

    if showLengths or showGaps or showNs or matchAmbiguous or not concise:
        result.append("""
<p>

Key to abbreviations:
  <ul>
    """)

        if showLengths:
            result.append('<li>L: sequence Length.</li>')

        if showGaps:
            result.append('<li>G: number of Gaps in sequence.</li>')

        if showNs:
            result.append('<li>N: number of N characters in sequence.</li>')

        if not concise:
            result.append('<li>IM: Identical nucleotide Matches.</li>')

        if matchAmbiguous:
            result.append('<li>AM: Ambiguous nucleotide Matches.</li>')

        result.append("""
    <li>GG: Gap/Gap matches (both sequences have gaps).</li>
    <li>G?: Gap/Non-gap mismatches (one sequence has a gap).</li>
    <li>NE: Non-equal nucleotide mismatches.</li>
  </ul>
</p>
""")

    return '\n'.join(result)


def dataCell(id1: str, id2: str, square: bool, readNumbers: Dict[str, int],
             upperOnly: bool) -> bool:
    """
    Should a cell have a value computed for it, or is it empty?

    @param id1: The row read id.
    @param id2: The column read id.
    @param square: If C{True} we are making a square table of a set of
        sequences against themselves (in which case we show nothing on the
        diagonal).
    @param readNumbers: A C{dict} mapping read ids to row numbers (only
        used if square is C{True} (in which case reads1 is the same as reads2).
    @param upperOnly: If C{True}, only compute values for the upper diagonal.
    @return: C{True} if the cell should be filled in, else C{False}.
    """
    if id1 == id2 and square:
        return False

    if not upperOnly:
        return True

    return readNumbers[id1] < readNumbers[id2]


def collectData(reads1, reads2, square, matchAmbiguous, pairwiseAlign,
                verbose, upperOnly=False):
    """
    Get pairwise matching statistics for two sets of reads.

    @param reads1: An C{OrderedDict} of C{str} read ids whose values are
        C{Read} instances. These will be the rows of the table.
    @param reads2: An C{OrderedDict} of C{str} read ids whose values are
        C{Read} instances. These will be the columns of the table.
    @param square: If C{True} we are making a square table of a set of
        sequences against themselves (in which case we show nothing on the
        diagonal).
    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct. Otherwise, we are strict
        and insist that only non-ambiguous nucleotides can contribute to the
        matching nucleotide count.
    @param pairwiseAlign: If C{True}, pairwise-align the sequences.
    @param verbose: If C{True}, print progress output.
    @param upperOnly: If C{True}, only compute values for the upper diagonal.
    """
    readNumbers = {}
    comparisons = 0
    for readNumber, id1 in enumerate(reads1):
        readNumbers[id1] = readNumber

    for id1 in reads1:
        for id2 in reads2:
            if dataCell(id1, id2, square, readNumbers, upperOnly):
                comparisons += 1

    result = defaultdict(dict)
    count = 0

    for id1, read1 in reads1.items():
        for id2, read2 in reads2.items():
            if dataCell(id1, id2, square, readNumbers, upperOnly):
                count += 1
                if pairwiseAlign:
                    r1, r2 = align([read1, read2], args)
                else:
                    r1, r2 = read1, read2
                if verbose:
                    print(f"Comparing {count}/{comparisons} {id1!r} "
                          f"and {id2!r}.", file=sys.stderr)
                match = compareDNAReads(r1, r2, matchAmbiguous=matchAmbiguous)
                # print(match, file=sys.stderr)
                if not matchAmbiguous:
                    assert match['match']['ambiguousMatchCount'] == 0
                result[id1][id2] = result[id2][id1] = match

    return result, readNumbers


def computeIdentity(read1, read2, stats, matchAmbiguous, read1Len, digits):
    """
    Compute nucleotide identity for two reads (as a fraction of the number
    of relevant nucleotides in the first read).

    @param read1: A C{Read} instance.
    @param read2: A C{Read} instance.
    @param stats: A C{dict} as returned by C{compareDNAReads}.
    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct. Otherwise, we are strict
        and insist that only non-ambiguous nucleotides can contribute to the
        matching nucleotide count.
    @param read1Len: The C{int} length of read1 (excluding gap characters).
    @param digits: The C{int} number of digits to round values to.
    """
    match = stats['match']
    if matchAmbiguous:
        numerator = match['identicalMatchCount'] + match['ambiguousMatchCount']
        denominator = read1Len
    else:
        numerator = match['identicalMatchCount']
        denominator = read1Len - len(stats['read1']['ambiguousOffsets'])

    # Note that the strict identity may be higher or lower than the
    # ambiguous identity even though an ambiguous match sounds like it
    # should always be more lenient and therefore result in a better
    # match. While it is true that the raw number of matched characters
    # will always increase when strict=False, the fraction of matched
    # characters may go down.
    #
    # The strict value can be higher because read1 might have many
    # ambiguous characters but very few of them may match read2. In that
    # case the overall fraction of matching characters will be pulled down
    # from the strict fraction when the ambiguous are included.
    #
    # Similarly, read1 may have many ambiguous characters, all of which are
    # matched by read2 and this can pull the overall identity higher than
    # the strict identity.
    #
    # If you want to play with this, change the following to True and
    # details of the case in which strictIdentity > ambiguousIdentity will
    # be printed to standard error.
    if False:
        # Apologies for the code duplication here. Above I only want to
        # compute one of the fractions, not both.
        numerator = match['identicalMatchCount'] + match['ambiguousMatchCount']
        denominator = read1Len
        ambiguousIdentity = numerator / denominator

        numerator = match['identicalMatchCount']
        denominator = read1Len - len(stats['read1']['ambiguousOffsets'])
        strictIdentity = numerator / denominator

        if strictIdentity > ambiguousIdentity:
            read1AmbiguityCount = len(stats['read1']['ambiguousOffsets'])
            print(f'Strict identity {strictIdentity} > ambiguous identity '
                  f'{ambiguousIdentity}. Read1 length {read1Len}. Read1 '
                  f'ambiguity count {read1AmbiguityCount}. Stats:\n{stats}',
                  file=sys.stderr)

    return round(numerator / denominator, digits)


def textTable(tableData, reads1, reads2, readNumbers, square, matchAmbiguous,
              gapChars, numberedColumns, upperOnly=False, digits=3,
              addZeroes=False):
    """
    Make a text table showing inter-sequence distances.

    @param tableData: A C{defaultdict(dict)} keyed by read ids, whose values
        are the dictionaries returned by compareDNAReads.
    @param reads1: An C{OrderedDict} of C{str} read ids whose values are
        C{Read} instances. These will be the rows of the table.
    @param reads2: An C{OrderedDict} of C{str} read ids whose values are
        C{Read} instances. These will be the columns of the table.
    @param readNumbers: A C{dict} mapping read ids to row numbers (only
        used if square is C{True} (in which case reads1 is the same as reads2).
    @param square: If C{True} we are making a square table of a set of
        sequences against themselves (in which case we show nothing on the
        diagonal).
    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct. Otherwise, we are strict
        and insist that only non-ambiguous nucleotides can contribute to the
        matching nucleotide count.
    @param gapChars: A C{str} of sequence characters considered to be gaps.
    @param numberedColumns: If C{True}, use (row) numbers for column names.
    @param upperOnly: If C{True}, only show values for the upper diagonal.
    @param digits: The C{int} number of digits to round identities to.b
    @param addZeroes: If C{True}, add trailing zeroes to identities so they
        all have the same width.
    """
    readLengths1 = getReadLengths(reads1.values(), gapChars)
    titles = ['ID']
    if numberedColumns:
        titles.extend(str(i + 1) for i in range(len(reads2)))

        if upperOnly and numberedColumns:
            titles.pop(1)
            titles[-1] = list(reads2)[-1]
    else:
        titles.extend(reads2)

    print('\t'.join(titles))

    for rowCount, (id1, read1) in enumerate(reads1.items(), start=1):
        if upperOnly and numberedColumns and rowCount == len(reads1):
            # We don't print the last row when only showing the upper
            # diagonal, because it will be empty. It's name will appear at
            # the top of the final column.
            continue
        prefix = f'{rowCount}: ' if numberedColumns else ''
        print(f'{prefix}{id1}', end='')
        for id2, read2 in reads2.items():
            if readNumbers[id2] == 0 and square:
                # The whole first column will be empty if we're making a
                # square array.
                continue
            if dataCell(id1, id2, square, readNumbers, upperOnly):
                identity = computeIdentity(
                    read1, read2, tableData[id1][id2], matchAmbiguous,
                    readLengths1[id1], digits)

                if addZeroes:
                    print(f'\t{identity:.{digits}f}', end='')
                else:
                    print(f'\t{identity}', end='')
            else:
                print('\t', end='')
        print()


def htmlTable(tableData, reads1, reads2, square, readNumbers, matchAmbiguous,
              colors, concise=False, showLengths=False, showGaps=False,
              showNs=False, footer=False, div=False, gapChars='-',
              numberedColumns=False, upperOnly=False, digits=3,
              addZeroes=False, highlightBest=False):
    """
    Make an HTML table showing inter-sequence distances.

    @param tableData: A C{defaultdict(dict)} keyed by read ids, whose values
        are the dictionaries returned by compareDNAReads.
    @param reads1: An C{OrderedDict} of C{str} read ids whose values are
        C{Read} instances. These will be the rows of the table.
    @param reads2: An C{OrderedDict} of C{str} read ids whose values are
        C{Read} instances. These will be the columns of the table.
    @param readNumbers: A C{dict} mapping read ids to row numbers (only
        used if square is C{True} (in which case reads1 is the same as reads2).
    @param square: If C{True} we are making a square table of a set of
        sequences against themselves (in which case we show nothing on the
        diagonal).
    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct. Otherwise, we are strict
        and insist that only non-ambiguous nucleotides can contribute to the
        matching nucleotide count.
    @param colors: A C{list} of (threshold, color) tuples, where threshold is a
        C{float} and color is a C{str} to be used as a cell background. This
        is as returned by C{parseColors}.
    @param concise: If C{True}, do not show match details.
    @param showLengths: If C{True}, include the lengths of sequences.
    @param showGaps: If C{True}, include the number of gaps in sequences.
    @param showGaps: If C{True}, include the number of N characters in
        sequences.
    @param footer: If C{True}, incude a footer row giving the same information
        as found in the table header.
    @param div: If C{True}, return an HTML <div> fragment only, not a full HTML
        document.
    @param gapChars: A C{str} of sequence characters considered to be gaps.
    @param numberedColumns: If C{True}, use (row) numbers for column names.
    @param upperOnly: If C{True}, only show values for the upper diagonal.
    @param digits: The C{int} number of digits to round identities to.b
    @param addZeroes: If C{True}, add trailing zeroes to identities so they
        all have the same width.
    @param highlightBest: If C{True}, highlight the best identity value
        in each row.
    @return: An HTML C{str} showing inter-sequence distances.
    """
    readLengths1 = getReadLengths(reads1.values(), gapChars)
    readLengths2 = getReadLengths(reads2.values(), gapChars)
    result = []
    append = result.append

    def writeHeader():
        # The header row of the table.
        append('    <tr>')
        append('    <td>&nbsp;</td>')
        for count, read2 in enumerate(reads2.values(), start=1):
            if count == 1 and square:
                # The first column will be empty, so skip it.
                continue
            append('    <td class="title"><span class="name">%s</span>' %
                   (count if (
                       upperOnly and numberedColumns and count != len(reads2))
                    else read2.id))
            if not square:
                if showLengths:
                    append('    <br>L:%d' % readLengths2[read2.id])
                if showGaps:
                    append('    <br>G:%d' % (
                        len(read2) - readLengths2[read2.id]))
                if showNs:
                    append('    <br>N:%d' % read2.sequence.count('N'))
            append('    </td>')
        append('    </tr>')

    if div:
        append('<div>')
    else:
        append('<!DOCTYPE HTML>')
        append('<html>')
        append('<head>')
        append('<meta charset="UTF-8">')
        append('</head>')
        append('<body>')

    append('<style>')
    append("""
        table {
            border-collapse: collapse;
        }
        table, td {
            border: 1px solid #ccc;
        }
        tr:hover {
            background-color: #f2f2f2;
        }
        td {
            vertical-align: top;
            font-size: 14px;
        }
        span.name {
            font-weight: bold;
        }
        span.best {
            font-weight: bold;
        }
        td.nt-identity {
            text-align: right;
        }
    """)

    # Add color style information for the identity thresholds.
    for threshold, color in colors:
        append('.%s { background-color: %s; }' % (
            thresholdToCssName(threshold), color))

    append('</style>')

    if not div:
        append(explanation(
            matchAmbiguous, concise, showLengths, showGaps, showNs))
    append('<div style="overflow-x:auto;">')
    append('<table>')
    append('  <tbody>')

    # Pre-process to find the best identities in each sample row.
    bestIdentityForId = {}
    identities = defaultdict(dict)

    for id1, read1 in reads1.items():
        # Look for best identity for the sample.
        bestIdentity = -1.0
        for id2, read2 in reads2.items():
            if dataCell(id1, id2, square, readNumbers, upperOnly):
                identity = computeIdentity(
                    read1, read2, tableData[id1][id2], matchAmbiguous,
                    readLengths1[id1], digits)
                identities[id1][id2] = identity
                if identity > bestIdentity:
                    bestIdentity = identity

        bestIdentityForId[id1] = bestIdentity

    writeHeader()

    # The main body of the table.
    for rowCount, (id1, read1) in enumerate(reads1.items(), start=1):
        if upperOnly and numberedColumns and rowCount == len(reads1):
            # We don't print the last row when only showing the upper
            # diagonal, because it will be empty. It's name will appear at
            # the top of the final column.
            continue

        read1Len = readLengths1[id1]
        append('    <tr>')
        append('      <td class="title"><span class="name">%s%s</span>' % (
            f"{rowCount}: " if numberedColumns else "", id1))
        if showLengths:
            append('<br/>L:%d' % read1Len)
        if showGaps:
            append('<br/>G:%d' % (len(read1) - read1Len))
        if showNs:
            append('<br/>N:%d' % read1.sequence.count('N'))
        append('</td>')
        for id2, read2 in reads2.items():
            if readNumbers[id2] == 0 and square:
                # The whole first column will be empty if we're making a
                # square array.
                continue

            if not dataCell(id1, id2, square, readNumbers, upperOnly):
                append('<td>&nbsp;</td>')
                continue

            identity = identities[id1][id2]

            append('      <td class="nt-identity %s">' % thresholdToCssName(
                thresholdForIdentity(identity, colors)))

            # The maximum percent identity.
            if highlightBest and identity == bestIdentityForId[id1]:
                scoreStyle = ' class="best"'
            else:
                scoreStyle = ''

            if addZeroes:
                append(f'<span{scoreStyle}>{identity:.{digits}f}</span>')
            else:
                append(f'<span{scoreStyle}>{identity}</span>')

            if not concise:
                match = tableData[id1][id2]['match']
                append('<br/>IM:%d' % match['identicalMatchCount'])

                if matchAmbiguous:
                    append('<br/>AM:%d' % match['ambiguousMatchCount'])

                append(
                    '<br/>GG:%d'
                    '<br/>G?:%d'
                    '<br/>NE:%d' %
                    (match['gapGapMismatchCount'],
                     match['gapMismatchCount'],
                     match['nonGapMismatchCount']))
            append('      </td>')
        append('    </tr>')

    if footer:
        writeHeader()

    append('  </tbody>')
    append('</table>')
    append('</div>')

    if div:
        append('</div>')
    else:
        append('</body>')
        append('</html>')

    return '\n'.join(result)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Print a FASTA sequence identity table.')

    parser.add_argument(
        '--text', action='store_true',
        help='If specified, just print a simple text table.')

    parser.add_argument(
        '--strict', action='store_true',
        help='If given, do not allow ambiguous nucleotide symbols to match.')

    parser.add_argument(
        '--numberedColumns', action='store_true',
        help=('Use a sequence (row) number as the header of each column '
              'instead of the sequence id.'))

    parser.add_argument(
        '--concise', action='store_true',
        help='If given, do not show match details.')

    parser.add_argument(
        '--digits', default=3, type=int,
        help='The number of digits to round identities to.')

    parser.add_argument(
        '--addZeroes', action='store_true',
        help=('If given, make sure all identities have the same number of '
              'digits shown by adding zeroes when needed.'))

    parser.add_argument(
        '--showLengths', action='store_true',
        help='If given, show the lengths of sequences.')

    parser.add_argument(
        '--showGaps', action='store_true',
        help='If given, show the number of gaps in sequences.')

    parser.add_argument(
        '--showNs', action='store_true',
        help=('If given, show the number of fully ambiguous N characters in '
              'sequences.'))

    parser.add_argument(
        '--footer', action='store_true',
        help='If given, also show sequence ids at the bottom of the table.')

    parser.add_argument(
        '--div', action='store_true',
        help=('If given, print an HTML <div> fragment only, not a full HTML '
              'document (ignored if --text is used).'))

    parser.add_argument(
        '--fastaFile2', type=open, metavar='FILENAME',
        help=('The name of a second FASTA input file. If no second FASTA '
              'file is given, sequences in the first FASTA file will be '
              'compared with each other.'))

    parser.add_argument(
        '--gapChars', default='-', metavar='CHARS',
        help=('The sequence characters that should be considered to be gaps. '
              'These characters will be ignored in computing sequence lengths '
              'and identity fractions.'))

    parser.add_argument(
        '--defaultColor', default='white',
        help=('The (background) color for cells. This will be used for all '
              'cells that do not otherwise have a color due to use of '
              '--color. This option is ignored if --text is given.'))

    parser.add_argument(
        '--color', action='append',
        help=('Specify cell background coloring. This option must be given as '
              'a space separated "value color" pair. The value is an identity '
              'fraction in [0..1] and the color is any color '
              'specification that can be given to CSS. This argument can be '
              'repeated. E.g., --color "0.9 red" --color "0.75 rgb(23, 190, '
              '207)" --color "0.1 #CF3CF3". Cells will be colored using the '
              'color of the highest identity fraction they satisfy. The '
              'default is to color all cells with the --defaultColor color. '
              'This option is ignored if --text is given.'))

    parser.add_argument(
        '--align', action='store_true',
        help='Do pairwise alignment of all sequences.')

    parser.add_argument(
        '--upperOnly', action='store_true',
        help=('Only show the upper diagonal (only valid when a single FASTA '
              'input file is given, resulting in a square all-against-all '
              'identity table).'))

    parser.add_argument(
        '--aligner', default='edlib', choices=('edlib', 'mafft', 'needle'),
        help='The alignment algorithm to use.')

    parser.add_argument(
        '--alignerOptions',
        help=('Optional arguments to pass to the alignment algorithm. If the '
              'aligner is mafft, the default options are %r. If needle, "%s". '
              'Do not try to set the number of threads here - use the '
              '--threads argument instead. If you are using mafft, see %s '
              'for some possible option combinations.' %
              (MAFFT_DEFAULT_ARGS, NEEDLE_DEFAULT_ARGS, MAFFT_ALGORITHMS_URL)))

    parser.add_argument(
        '--verbose', action='store_true',
        help='Print progress to standard error.')

    parser.add_argument(
        '--highlightBest', action='store_true',
        help='Highlight the best identity in each row.')

    addFASTACommandLineOptions(parser)
    addFASTAFilteringCommandLineOptions(parser)
    addFASTAEditingCommandLineOptions(parser)

    args = parser.parse_args()

    colors = parseColors(args.color, args.defaultColor)
    # Sanity check - the last threshold must be zero.
    assert colors[-1][0] == 0.0

    reads = parseFASTAEditingCommandLineOptions(
        args, parseFASTAFilteringCommandLineOptions(
            args, parseFASTACommandLineOptions(args)))

    # Collect the reads into a dict, keeping the insertion order.
    reads1 = OrderedDict()
    for read in reads:
        reads1[read.id] = read

    if args.fastaFile2:
        if args.upperOnly:
            print('The --upperOnly option is not supported if give two FASTA '
                  'input files.', file=sys.stderr)
            sys.exit(1)

        square = False
        reads2 = OrderedDict()
        # The next line is a total hack, to trick parseFASTACommandLineOptions
        # into reading a second FASTA file.
        args.fastaFile = args.fastaFile2
        for read in parseFASTAFilteringCommandLineOptions(
                args, parseFASTACommandLineOptions(args)):
            reads2[read.id] = read
    else:
        square = True
        reads2 = reads1

    matchAmbiguous = not args.strict
    tableData, readNumbers = collectData(
        reads1, reads2, square, matchAmbiguous, args.align, args.verbose,
        args.upperOnly)

    if args.text:
        textTable(tableData, reads1, reads2, readNumbers, square,
                  matchAmbiguous, args.gapChars, args.numberedColumns,
                  upperOnly=args.upperOnly, digits=args.digits,
                  addZeroes=args.addZeroes)
    else:
        print(
            htmlTable(
                tableData, reads1, reads2, square, readNumbers, matchAmbiguous,
                colors=colors, concise=args.concise,
                showLengths=args.showLengths, showGaps=args.showGaps,
                showNs=args.showNs, footer=args.footer, div=args.div,
                gapChars=args.gapChars, numberedColumns=args.numberedColumns,
                upperOnly=args.upperOnly, digits=args.digits,
                addZeroes=args.addZeroes, highlightBest=args.highlightBest))
