import sys
from os.path import join
from subprocess import CalledProcessError
from tempfile import mkdtemp
from shutil import rmtree

import edlib

from dark.fasta import FastaReads
from dark.process import Executor
from dark.reads import Reads


def mafft(reads, verbose=False, options=None, threads=None):
    """
    Run a MAFFT alignment and return the sequences.

    @param reads: An iterable of multiple reads.
    @param verbose: If C{True} print progress info to sys.stderr.
    @param options: A C{str} of options to pass to mafft.
    @return: A C{Reads} instance with the aligned sequences.
    """
    tempdir = mkdtemp()

    infile = join(tempdir, 'input.fasta')
    out = join(tempdir, 'result.fasta')

    Reads(reads).save(infile)

    if verbose:
        print('Running mafft.', file=sys.stderr)

    Executor().execute("mafft %s %s '%s' > '%s'" % (
        ('' if threads is None else '--thread %d' % threads),
        options or '', infile, out))

    # Use 'list' in the following to force reading the FASTA from disk.
    result = Reads(list(FastaReads(out)))
    rmtree(tempdir)

    return result


def needle(reads, verbose=False, options=None):
    """
    Run a Needleman-Wunsch alignment and return the two sequences.

    @param reads: An iterable of two reads.
    @param verbose: If C{True} print progress info to sys.stderr.
    @param options: Additional options to pass to needle.
    @return: A C{Reads} instance with the two aligned sequences.
    """
    tempdir = mkdtemp()

    file1 = join(tempdir, 'file1.fasta')
    with open(file1, 'w') as fp:
        print(reads[0].toString('fasta'), end='', file=fp)

    file2 = join(tempdir, 'file2.fasta')
    with open(file2, 'w') as fp:
        print(reads[1].toString('fasta'), end='', file=fp)

    out = join(tempdir, 'result.fasta')

    def useStderr(e):
        return "Sequences too big. Try 'stretcher'" not in e.stderr

    if verbose:
        print('Running needle.', file=sys.stderr)
    try:
        Executor().execute(
            "needle -asequence '%s' -bsequence '%s' %s "
            "-outfile '%s' -aformat fasta" % (
                file1, file2, options or '', out), useStderr=useStderr)
    except CalledProcessError as e:
        if useStderr(e):
            raise
        else:
            if verbose:
                print('Sequences too long for needle. Falling back to '
                      'stretcher. Be patient!', file=sys.stderr)
            Executor().execute("stretcher -asequence '%s' -bsequence '%s' "
                               "-auto "
                               "-outfile '%s' -aformat fasta" % (
                                   file1, file2, out))

    # Use 'list' in the following to force reading the FASTA from disk.
    result = Reads(list(FastaReads(out)))
    rmtree(tempdir)

    return result


def edlibAlign(reads, gapSymbol='-', strict=True):
    """
    Run an edlib alignment and return the sequences.

    @param reads: An iterable of at least two reads.
    @param gapSymbol: A C{str} 1-character symbol to use for gaps.
    @param strict: Ensure that C{reads} only has two reads.
    @raise ValueError: If C{strict} is C{True} and there are more than two
       reads passed.
    @return: A C{Reads} instance with the aligned sequences.
    """
    r1, r2, *rest = list(reads)

    if strict and len(rest):
        raise ValueError(f'{len(rest)} unexpected extra arguments')

    alignment = edlib.getNiceAlignment(
        edlib.align(r1.sequence, r2.sequence, mode='NW', task='path'),
        r1.sequence, r2.sequence, gapSymbol=gapSymbol)

    return Reads([
        r1.__class__(r1.id, alignment['query_aligned']),
        r2.__class__(r2.id, alignment['target_aligned']),
    ])
