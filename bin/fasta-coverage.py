#!/usr/bin/env python

from collections import defaultdict

from dark.dna import AMBIGUOUS
from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Print coverage information to standard output.')

    parser.add_argument(
        '--chars', '-c', default='N',
        help='The characters that count as non-covered positions.')

    parser.add_argument(
        '--noNames', dest='printNames', action='store_false',
        help='Do not print the ids of input sequences.')

    parser.add_argument(
        '--strict', action='store_true',
        help=('Treat all ambiguous nucleotides as not covered. Additional '
              'characters specified with --chars will also be counted as not '
              'covered.'))

    parser.add_argument(
        '--gc', action='store_true',
        help='Include the GC percentage in the output.')

    parser.add_argument(
        '--sep', default='\t',
        help='The output separator character.')

    parser.add_argument(
        '--printCharCounts', action='store_true',
        help='Print the count for each character found in each read.')

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    notCoveredCount = 0
    notCovered = set(args.chars)

    if args.strict:
        notCovered.update(
            code for code, chars in AMBIGUOUS.items() if len(chars) > 1)

    for read in reads:
        counts = defaultdict(int)
        for char in read.sequence:
            notCoveredCount += (char in notCovered)
            counts[char] += 1

        line = [read.id] if args.printNames else []

        if args.printCharCounts:
            line.extend(f'{char}:{count}'
                        for char, count in sorted(counts.items()))

        length = len(read)

        if args.gc:
            gc = 0
            for char in 'GC':
                if char in counts:
                    gc += counts[char]

            line.append(f'gc:{gc / length:.3f}')

        coverage = (length - notCoveredCount) / length if length else 0.0

        line.extend((
            f'length:{length}',
            f'not-covered:{notCoveredCount}',
            f'coverage:{coverage:.3f}',
        ))

        print(args.sep.join(line))
