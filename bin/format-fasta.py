#!/usr/bin/env python

from __future__ import print_function

import sys
import hashlib

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given FASTA on stdin, write it to stdout according to a '
                     'given format.'))

    parser.add_argument(
        '--prefix',
        help=('The string prefix to be used. E.g., use --prefix r for raw '
              'strings. The default is "f" (i.e., use Python f-strings) for '
              'Python 3.6 and later, or else "" (i.e., regular strings).'))

    parser.add_argument(
        '--start', default=1, type=int,
        help='The integer to start counting from, for the "count" variable.')

    parser.add_argument(
        '--end',
        help=('A string to use as the "end" keyword argument for print(). '
              'This will be passed to eval, allowing convenient use of '
              'Python string constructs such as \\t on the command line. '
              'If not specified, a newline will be printed.'))

    parser.add_argument(
        '--format',
        help=("The output format. The default is to reproduce the input. "
              "The format can be that used by Python to specify items in a "
              'dictionary. Available dictionary items are "id", "sequence", '
              '"quality", "length", "count", and "md5" (of the sequence). A '
              "dark.reads.Read instance (named 'read') and an 'md5' "
              "function (that takes and returns a string) are both also in "
              "scope for each FASTA record and these can be accessed using "
              "f-string format (with Python 3.6 or greater). The passed "
              "format string will be eval'd in a double quoted string, as "
              "with the --end string, allowing the use of Python string "
              "conveniences such as \\t on the command line."))

    addFASTACommandLineOptions(parser)
    args = parser.parse_args()
    reads = parseFASTACommandLineOptions(args)

    format_ = args.format
    end = '\n' if args.end is None else eval('"' + args.end + '"')

    if args.prefix is None:
        prefix = '' if sys.version_info < (3, 6) else 'f'
    else:
        prefix = args.prefix

    if format_ is None:
        format_ = (r'>%(id)s\n%(sequence)s' if args.fasta else
                   r'@%(id)s\n%(sequence)s\n+\n%(quality)s')

    needMd5 = '%(md5)' in format_
    needLength = '%(length)' in format_

    def md5(s):
        return hashlib.md5(s.encode('utf-8')).hexdigest()

    for count, read in enumerate(reads, start=args.start):
        sequence = read.sequence
        d = {
            'count': count,
            'id': read.id,
            'quality': read.quality,
            'sequence': sequence,
        }
        if needMd5:
            d['md5'] = md5(sequence)
        if needLength:
            d['length'] = len(sequence)

        print(eval(prefix + '"' + format_ + '" % d'), end=end)
