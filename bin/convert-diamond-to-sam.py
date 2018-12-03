#!/usr/bin/env python

# See https://samtools.github.io/hts-specs/SAMv1.pdf for the SAM file
# format specification.

from __future__ import print_function, division

import sys
import argparse
from os.path import basename
from tempfile import TemporaryFile
from functools import partial

from dark import __version__ as VERSION
from dark.btop import btop2cigar
from dark.diamond.conversion import diamondTabularFormatToDicts, FIELDS
from dark.reads import DNARead

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Convert DIAMOND tabular format to SAM. The DIAMOND '
                 'invocation *must* include --outfmt 6 %s' % FIELDS))

parser.add_argument(
    '--printFields', default=False, action='store_true',
    help=('Print the field names in the order that they must be given to '
          'diamond --outfmt 6 to produce correct input for this script, '
          'then exit.'))

parser.add_argument(
    '--mappingQuality', type=int, default=255,
    help=('The mapping quality to use for MAPQ (field 5). The default (255) '
          'indicates that mapping quality information is not available.'))

parser.add_argument(
    '--ram', action='store_true', default=False,
    help=('Do not use a temporary file to hold the non-header SAM output. '
          'This will run faster but use more memory since all non-header SAM '
          'output will be stored in RAM and only written out when the full '
          'header can be determined.'))

parser.add_argument(
    '--keepDescriptions', action='store_true', default=False,
    help=('Do not discard text after the first space in query or subject '
          'sequence ids. Note that. Note that this violates the SAM '
          'specification, but since SAM files are TAB-separated there '
          'is probably only a small chance this will cause any problems '
          'downstream.'))

args = parser.parse_args()

if args.printFields:
    print(FIELDS)
    sys.exit(0)

idOnly = not args.keepDescriptions
mappingQuality = args.mappingQuality
ram = args.ram

if 0 > mappingQuality > 255:
    raise ValueError('Mapping quality must be between 0 and 255 (inclusive)')

referenceLengths = {}

if ram:
    nonHeaderLines = []
    emit = nonHeaderLines.append
else:
    tf = TemporaryFile(mode='w+t', encoding='utf-8')
    emit = partial(print, file=tf)

for match in diamondTabularFormatToDicts(sys.stdin):
    qseqid = match['qseqid'].split()[0] if idOnly else match['qseqid']
    stitle = match['stitle'].split()[0] if idOnly else match['stitle']

    # The subject length is ALWAYS in amino acids in DIAMOND.
    referenceLengths[stitle] = 3 * match['slen']

    # If the query frame is less than zero, the match was with a reverse
    # complemented translation of the query. Put the reverse compliment
    # into the SAM output, which seems to be standard / accepted practice
    # based on my web searches. See e.g., https://www.biostars.org/p/131891/
    # for what Bowtie2 does and for some comments on this issue for SAM/BAM
    # files in general.
    if match['qframe'] > 0:
        flag = 0
        qseq = match['qseq']
        qqual = match['qqual'] or '*'
    else:
        flag = 16
        qseq = DNARead('id', match['qseq']).reverseComplement().sequence
        qqual = match['qqual'][::-1] if match['qqual'] else '*'

    # Make a CIGAR string, including hard-clipped bases at the start and
    # end of the query (DIAMOND outputs a hard-clipped query sequence).
    startClipCount = match['qstart'] - 1
    endClipCount = match['qlen'] - match['qend']

    assert startClipCount >= 0
    assert endClipCount >= 0, (
        'Query sequence %s has length %d but the qend value is %d' %
        (qseq, len(match['qseq']), match['qend']))

    cigar = (
        ('%dH' % startClipCount if startClipCount else '') +
        btop2cigar(match['btop'], concise=False, aa=True) +
        ('%dH' % endClipCount if endClipCount else ''))

    emit('\t'.join(map(str, [
        # 1. QNAME
        qseqid,
        # 2. FLAG
        flag,
        # 3. RNAME
        stitle,
        # 4. POS. This needs to be a 1-based offset into the
        # nucleotide-equivalent of the DIAMOND subject sequence (which was
        # a protein since that is how DIAMOND operates). Because DIAMOND
        # gives back a 1-based protein location, we adjust to 0-based,
        # multiply by 3 to get to nucleotides, then adjust to 1-based.
        3 * (match['sstart'] - 1) + 1,
        # 5. MAPQ
        mappingQuality,
        # 6. CIGAR
        cigar,
        # 7. RNEXT
        '*',
        # 8. PNEXT
        0,
        # 9. TLEN
        0,
        # 10. SEQ
        qseq,
        # 11. QUAL
        qqual,
        # 12. Alignment score
        'AS:i:%d' % int(match['bitscore'])])))


progName = basename(sys.argv[0])

# Print SAM headers.
print('\n'.join(
    [
        '@PG\tID:DIAMOND\tPN:DIAMOND',
        '@PG\tID:%s\tPN:%s (version %s)\tCL:%s %s\tPP:DIAMOND' %
        (progName, progName, VERSION, progName, ' '.join(sys.argv[1:])),
        '@CO\t%s is from the dark-matter package '
        '(https://github.com/acorg/dark-matter/)' % progName,
    ] +
    [
        '@SQ\tSN:%s\tLN:%d' % (name, referenceLengths[name])
        for name in sorted(referenceLengths)
    ]))

# Print non-header lines.
if ram:
    print('\n'.join(nonHeaderLines))
else:
    tf.seek(0)
    for line in tf:
        print(line, end='')
    tf.close()
