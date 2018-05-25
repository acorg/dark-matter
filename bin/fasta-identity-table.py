#!/usr/bin/env python

from __future__ import print_function, division

import argparse
from collections import OrderedDict, defaultdict

from dark.dna import compareDNAReads
from dark.filter import (
    addFASTAFilteringCommandLineOptions, parseFASTAFilteringCommandLineOptions)
from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


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


def explanation(matchAmbiguous, concise, showLengths, showGaps):
    """
    Make an explanation of the output HTML table.

    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct. Otherwise, we are strict
        and insist that only non-ambiguous nucleotides can contribute to the
        matching nucleotide count.
    @param concise: If C{True}, do not show match detail abbreviations.
    @param showLengths: If C{True}, include the lengths of sequences.
    @param showGaps: If C{True}, include the number of gaps in sequences.
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

    if showLengths or showGaps or matchAmbiguous or not concise:
        result.append("""
<p>

Key to abbreviations:
  <ul>
    """)

        if showLengths:
            result.append('<li>L: sequence Length.</li>')

        if showGaps:
            result.append('<li>G: number of Gaps in sequence.</li>')

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


def collectData(reads1, reads2, square, matchAmbiguous):
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
    """
    result = defaultdict(dict)
    for id1, read1 in reads1.items():
        for id2, read2 in reads2.items():
            if id1 != id2 or not square:
                match = compareDNAReads(
                    read1, read2, matchAmbiguous=matchAmbiguous)['match']
                if not matchAmbiguous:
                    assert match['ambiguousMatchCount'] == 0
                result[id1][id2] = result[id2][id1] = match

    return result


def simpleTable(tableData, reads1, reads2, square, matchAmbiguous, gapChars):
    """
    Make a text table showing inter-sequence distances.

    @param tableData: A C{defaultdict(dict)} keyed by read ids, whose values
        are the dictionaries returned by compareDNAReads.
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
    @param gapChars: A C{str} of sequence characters considered to be gaps.
    """
    readLengths1 = getReadLengths(reads1.values(), gapChars)
    print('ID\t' + '\t'.join(reads2))

    for id1, read1 in reads1.items():
        read1Len = readLengths1[id1]
        print(id1, end='')
        for id2, read2 in reads2.items():
            if id1 == id2 and square:
                print('\t', end='')
            else:
                stats = tableData[id1][id2]
                identity = (
                    stats['identicalMatchCount'] +
                    (stats['ambiguousMatchCount'] if matchAmbiguous else 0)
                ) / read1Len
                print('\t%.4f' % identity, end='')
        print()


def htmlTable(tableData, reads1, reads2, square, matchAmbiguous, concise=False,
              showLengths=False, showGaps=False, footer=False, div=False,
              gapChars='-'):
    """
    Make an HTML table showing inter-sequence distances.

    @param tableData: A C{defaultdict(dict)} keyed by read ids, whose values
        are the dictionaries returned by compareDNAReads.
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
    @param concise: If C{True}, do not show match details.
    @param showLengths: If C{True}, include the lengths of sequences.
    @param showGaps: If C{True}, include the number of gaps in sequences.
    @param footer: If C{True}, incude a footer row giving the same information
        as found in the table header.
    @param div: If C{True}, return an HTML <div> fragment only, not a full HTML
        document.
    @param gapChars: A C{str} of sequence characters considered to be gaps.
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
        for read2 in reads2.values():
            append('    <td class="title"><span class="name">%s</span>' %
                   read2.id)
            if showLengths and not square:
                append('    <br>L:%d' % readLengths2[read2.id])
            if showGaps and not square:
                append('    <br>G:%d' % (len(read2) - readLengths2[read2.id]))
            append('    </td>')
        append('    </tr>')

    if div:
        append('<div>')
    else:
        append('<html>')
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
    """)
    append('</style>')

    if not div:
        append(explanation(matchAmbiguous, concise, showLengths, showGaps))
    append('<div style="overflow-x:auto;">')
    append('<table>')
    append('  <tbody>')

    # Pre-process to find the best identities in each sample row.
    bestIdentityForId = {}

    for id1, read1 in reads1.items():
        # Look for best identity for the sample.
        read1Len = readLengths1[id1]
        bestIdentity = -1.0
        for id2, read2 in reads2.items():
            if id1 != id2 or not square:
                stats = tableData[id1][id2]
                identity = (
                    stats['identicalMatchCount'] +
                    (stats['ambiguousMatchCount'] if matchAmbiguous else 0)
                ) / read1Len
                if identity > bestIdentity:
                    bestIdentity = identity
        bestIdentityForId[id1] = bestIdentity

    writeHeader()

    # The main body of the table.
    for id1, read1 in reads1.items():
        read1Len = readLengths1[id1]
        append('    <tr>')
        append('      <td class="title"><span class="name">%s</span>' % id1)
        if showLengths:
            append('<br/>L:%d' % read1Len)
        if showGaps:
            append('<br/>G:%d' % (len(read1) - read1Len))
        append('</td>')
        for id2, read2 in reads2.items():
            if id1 == id2 and square:
                append('<td>&nbsp;</td>')
                continue

            stats = tableData[id1][id2]
            identity = (
                stats['identicalMatchCount'] +
                (stats['ambiguousMatchCount'] if matchAmbiguous else 0) +
                stats['gapGapMismatchCount']
            ) / read1Len

            append('      <td>')

            # The maximum percent identity.
            if identity == bestIdentityForId[id1]:
                scoreStyle = ' class="best"'
            else:
                scoreStyle = ''
            append('<span%s>%.4f</span>' % (scoreStyle, identity))

            if not concise:
                append('<br/>IM:%d' % stats['identicalMatchCount'])

                if matchAmbiguous:
                    append('<br/>AM:%d' % stats['ambiguousMatchCount'])

                append(
                    '<br/>GG:%d'
                    '<br/>G?:%d'
                    '<br/>NE:%d' %
                    (stats['gapGapMismatchCount'],
                     stats['gapMismatchCount'],
                     stats['nonGapMismatchCount']))
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
        '--text', default=False, action='store_true',
        help='If specified, just print a simple text table')

    parser.add_argument(
        '--strict', default=False, action='store_true',
        help='If given, do not allow ambiguous nucleotide symbols to match')

    parser.add_argument(
        '--concise', default=False, action='store_true',
        help='If given, do not show match details')

    parser.add_argument(
        '--showLengths', default=False, action='store_true',
        help='If given, show the lengths of sequences')

    parser.add_argument(
        '--showGaps', default=False, action='store_true',
        help='If given, show the number of gaps in sequences')

    parser.add_argument(
        '--footer', default=False, action='store_true',
        help='If given, also show sequence ids at the bottom of the table')

    parser.add_argument(
        '--div', default=False, action='store_true',
        help=('If given, print an HTML <div> fragment only, not a full HTML '
              'document (ignored if --text is used)'))

    parser.add_argument(
        '--fastaFile2', type=open, metavar='FILENAME',
        help=('The name of a second FASTA input file. If no second FASTA '
              'file is given, sequences in the first FASTA file will be '
              'compared with each other.'))

    parser.add_argument(
        '--gapChars', default='-', metavar='CHARS',
        help=('The sequence characters that should be considered to be gaps. '
              'These characters will be ignored in computing sequence lengths '
              'and identity fractions'))

    addFASTACommandLineOptions(parser)
    addFASTAFilteringCommandLineOptions(parser)
    args = parser.parse_args()

    # Collect the reads into a dict, keeping the insertion order.
    reads1 = OrderedDict()
    for read in parseFASTAFilteringCommandLineOptions(
            args, parseFASTACommandLineOptions(args)):
        reads1[read.id] = read

    if args.fastaFile2:
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
    tableData = collectData(reads1, reads2, square, matchAmbiguous)

    if args.text:
        simpleTable(tableData, reads1, reads2, square, matchAmbiguous)
    else:
        print(
            htmlTable(tableData, reads1, reads2, square, matchAmbiguous,
                      concise=args.concise, showLengths=args.showLengths,
                      showGaps=args.showGaps, footer=args.footer, div=args.div,
                      gapChars=args.gapChars))
