#!/usr/bin/env python

from __future__ import print_function, division

import argparse
from collections import OrderedDict, defaultdict

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions
from dark.dna import compareDNAReads


def explanation(matchAmbiguous, concise, showLengths):
    """
    Make an explanation of the output HTML table.

    @param matchAmbiguous: If C{True}, count ambiguous nucleotides that are
        possibly correct as actually being correct. Otherwise, we are strict
        and insist that only non-ambiguous nucleotides can contribute to the
        matching nucleotide count.
    @param concise: If C{True}, do not show match detail abbreviations.
    @param showLengths: If C{True}, include the lengths of sequences.
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
        result.append('nucleotides, plus the number of ambiguously matching '
                      'nucleotides, ')
    else:
        result.append('nucleotides')

    result.append("""
plus the number of sites where both sequences have gaps. The denominator
is the length of the sequence <em>for the row</em>.

</p>
    """)

    if showLengths or matchAmbiguous or not concise:
        result.append("""
<p>

Key to abbreviations:
  <ul>
    """)

        if showLengths:
            result.append('<li>L: sequence Length.</li>')

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


def simpleTable(tableData, reads1, reads2, square, matchAmbiguous):
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
    """
    print('ID\t' + '\t'.join(reads2))

    for id1, read1 in reads1.items():
        read1Len = len(read1)
        print(id1, end='')
        for id2, read2 in reads2.items():
            if id1 == id2 and square:
                print('\t', end='')
            else:
                stats = tableData[id1][id2]
                identity = (
                    stats['identicalMatchCount'] +
                    (stats['ambiguousMatchCount'] if matchAmbiguous else 0) +
                    stats['gapGapMismatchCount']
                ) / read1Len
                print('\t%.4f' % identity, end='')
        print()


def htmlTable(tableData, reads1, reads2, square, matchAmbiguous, concise=False,
              showLengths=False, footer=False, div=False):
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
    @param footer: If C{True}, incude a footer row giving the same information
        as found in the table header.
    @param div: If C{True}, return an HTML <div> fragment only, not a full HTML
        document.
    @return: An HTML C{str} showing inter-sequence distances.
    """
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
                append('    <br>L:%d' % len(read2))
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
        append(explanation(matchAmbiguous, concise, showLengths))
    append('<div style="overflow-x:auto;">')
    append('<table>')
    append('  <tbody>')

    # Pre-process to find the best identities in each sample row.
    bestIdentityForId = {}

    for id1, read1 in reads1.items():
        # Look for best identity for the sample.
        read1Len = len(read1)
        bestIdentity = -1.0
        for id2, read2 in reads2.items():
            if id1 != id2 or not square:
                stats = tableData[id1][id2]
                identity = (
                    stats['identicalMatchCount'] +
                    (stats['ambiguousMatchCount'] if matchAmbiguous else 0) +
                    stats['gapGapMismatchCount']
                ) / read1Len
                if identity > bestIdentity:
                    bestIdentity = identity
        bestIdentityForId[id1] = bestIdentity

    writeHeader()

    # The main body of the table.
    for id1, read1 in reads1.items():
        read1Len = len(read1)
        append('    <tr>')
        append('      <td class="title"><span class="name">%s</span>' % id1)
        if showLengths:
            append('<br/>L:%d' % read1Len)
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

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()

    # Collect the reads into a dict, keeping the insertion order.
    reads1 = OrderedDict()
    for read in parseFASTACommandLineOptions(args):
        reads1[read.id] = read

    if args.fastaFile2:
        square = False
        reads2 = OrderedDict()
        # This is a total hack, to trick parseFASTACommandLineOptions into
        # also reading the second FASTA file.
        args.fastaFile = args.fastaFile2
        for read in parseFASTACommandLineOptions(args):
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
                      footer=args.footer, div=args.div))
