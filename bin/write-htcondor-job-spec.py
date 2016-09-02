#!/usr/bin/env python

"""
See the 'EPILOG' variable below, or (better) run with --help for help.
"""

from __future__ import print_function

import os
import sys
from Bio import SeqIO

DEFAULT_BLAST_ARGS = ''
DEFAULT_BLAST_DB = 'nt'
DEFAULT_BLAST_DB_DIR = '/usr/local/dark-matter/blast-dbs'
DEFAULT_EMAIL = 'tcj25@cam.ac.uk'
DEFAULT_BLAST_EXECUTABLE_DIR = '/usr/local/dark-matter/blast/bin'
DEFAULT_BLAST_EXECUTABLE_NAME = 'blastn'
DEFAULT_SEQUENCES_PER_BLAST = 100

EPILOG = """Given a FASTA file argument, write out the following:

  1) Files named 0.fasta, 1.fasta, 2.fasta, etc. each containing a maximum
     number of sequences (given by --seqs-per-blast).

  2) An HTCondor job spec file 'job.htcondor' that can be given to
     condor_submit to have it run BLAST on all the small FASTA files.

  3) A 'redo.sh' script that can be used to re-submit sub-FASTA files on which
     the initial processing failed.

  4) A 'post-process.sh' script that can be used to convert BLAST XML to our
     JSON format in the case that the last step of processing fails.

  5) A 'finalize.sh' script that will find all zero length .json result
     files and post-process them (if they have an XML file) and that also
     removes any empty .error files.

NOTE: files are overwritten. It is suggested you run this command in an
empty directory.
"""


def splitFASTA(params):
    """
    Read the FASTA file named params['fastaFile'] and print out its
    sequences into files named 0.fasta, 1.fasta, etc. with
    params['seqsPerJob'] sequences per file.
    """
    assert params['fastaFile'][-1] == 'a', ('You must specify a file in '
                                            'fasta-format that ends in '
                                            '.fasta')

    fileCount = count = seqCount = 0
    outfp = None
    with open(params['fastaFile']) as infp:
        for seq in SeqIO.parse(infp, 'fasta'):
            seqCount += 1
            if count == params['seqsPerJob']:
                outfp.close()
                count = 0
            if count == 0:
                outfp = open('%d.fasta' % fileCount, 'w')
                fileCount += 1
            count += 1
            outfp.write('>%s\n%s\n' % (seq.description, str(seq.seq)))
    outfp.close()
    return fileCount, seqCount


def printJobSpec(params):
    """
    Write out a job spec file for HTCondor to process all the small
    FASTA input files via BLAST and our JSON post-processor.
    """
    with open('job.htcondor', 'w') as outfp:
        outfp.write("""\
universe                  = vanilla
executable                = process.sh
should_transfer_files     = YES
when_to_transfer_output   = ON_EXIT
notify_user               = %(email)s
max_transfer_input_mb     = -1
max_transfer_output_mb    = -1
log                       = job.log

# Job summary:
#   FASTA input from %(fastaFile)s
#   %(sequenceCount)d sequences split into %(nJobs)d jobs of \
%(seqsPerJob)d sequences each.

arguments                 = $(Process) %(db)s %(blastArgs)s
input                     = $(Process).fasta
output                    = $(Process).done
error                     = $(Process).error
dont_encrypt_input_files  = $(Process).fasta
dont_encrypt_output_files = $(Process).json.bz2

queue %(nJobs)d
""" % params)


def printRedoScript(params):
    """
    Write out a shell script that writes a job spec file for HTCondor to
    process a single FASTA input file via BLAST and our JSON post-processor,
    runs that job, and removes the one-time spec file it wrote.
    """
    with open('redo.sh', 'w') as outfp:
        outfp.write("""\
#!/bin/sh -e

case $# in
    0) echo "Usage: `basename $0` jobid1, jobid2, ..." >&2; exit 1;;
esac

tmp=redo.tmp.$$
trap "rm -f $tmp" 0 1 2 3 15

    cat >$tmp <<EOF
universe                  = vanilla
executable                = process.sh
should_transfer_files     = YES
when_to_transfer_output   = ON_EXIT
notify_user               = %(email)s
max_transfer_input_mb     = -1
max_transfer_output_mb    = -1
log                       = job.log
EOF

for jobid in "$@"
do
    cat >>$tmp <<EOF

arguments                 = $jobid %(db)s %(blastArgs)s
input                     = $jobid.fasta
output                    = $jobid.done
error                     = $jobid.error
dont_encrypt_input_files  = $jobid.fasta
dont_encrypt_output_files = $jobid.json.bz2
queue
EOF

    rm -f $jobid.json.bz2 $jobid.error $jobid.done
done

rm -f job.log

condor_submit $tmp
""" % params)

    # Make the script executable so we can run it.
    os.chmod('redo.sh', 0o755)


def printProcessScript(params):
    """
    Write out a simple process script to call BLAST and convert its XML to
    our compressed JSON format.
    """
    with open('process.sh', 'w') as outfp:
        outfp.write("""\
#!/bin/sh

DM=/usr/local/dark-matter
export PYTHONPATH=$DM/dark-matter
export BLASTDB="%(dbDir)s"

jobid=$1
shift

db="$1"
shift

errs=$jobid.error

%(executableDir)s/%(executableName)s -query $jobid.fasta -db "$db" \
-outfmt 5 "$@" |
$DM/virtualenv/bin/python $DM/dark-matter/bin/convert-blast-xml-to-json.py \
--bzip2 > $jobid.json.bz2 2> $errs

if [ -s $errs ]
then
    echo "Completed WITH ERRORS ($errs) on `hostname` at `date`." > $jobid.done
else
    rm $errs
    echo "Completed on `hostname` at `date`." > $jobid.done
fi
""" % params)

    # Make the script executable so we can run it.
    os.chmod('process.sh', 0o755)


def printFinalizeScript(params):
    """
    Write out a script that post-processes any empty .json files that
    have non-empty .xml files, figures out which jobs should be re-run, and
    also removes zero-length error files.

    Note that we need bash in order to set the nullglob shell option. That
    prevents an error if there are no *.fasta files.
    """
    with open('finalize.sh', 'w') as outfp:
        outfp.write("""\
#!/usr/bin/env bash

shopt -s nullglob

# redo will hold the numbers of jobs that need re-running, if any.
redo=

for i in *.fasta
do
    n=`echo $i | cut -f1 -d.`
    error=$n.error
    json=$n.json.bz2

    if [ -f $error ]
    then
        if [ -s $error ]
        then
            echo "WARNING: $error is non-empty." >&2
        else
            rm -f $error
        fi
    fi

    if [ -f $json ]
    then
        size=`wc -c < $json | awk '{print $1}'`
        # bzip2 run on empty input results in a file of 14 bytes.
        if [ $size -eq 14 ]
        then
            echo "WARNING: $json appears to be empty." >&2
        fi
    else
        echo "WARNING: $json does not exist. Job $n should be re-run." >&2
        redo="$redo $n"
    fi
done

if [ -n "$redo" ]
then
    echo "Some jobs must be re-run. Please execute the following:"
    echo "./redo.sh $redo"
fi
""" % params)

    # Make the script executable so we can run it.
    os.chmod('finalize.sh', 0o755)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Given a FASTA file, write an HTCondor job spec for BLAST',
        epilog=EPILOG)
    parser.add_argument(
        'fasta', metavar='FASTA-file',
        help='the FASTA file of sequences to BLAST.')
    parser.add_argument(
        '--seqs-per-blast', metavar='N',
        type=int, default=DEFAULT_SEQUENCES_PER_BLAST, dest='seqsPerJob',
        help='the number (>0) of sequences to pass to BLAST on each run.')
    parser.add_argument(
        '--blast-db-name', metavar='BLAST-database-name',
        default=DEFAULT_BLAST_DB, dest='db',
        help='the BLAST database to run against.')
    parser.add_argument(
        '--email', metavar='name@host',
        default=DEFAULT_EMAIL, dest='email',
        help='the email address to send the job completed message to.')
    parser.add_argument(
        '--blast-executable-name', metavar='BLAST-executable-name',
        default=DEFAULT_BLAST_EXECUTABLE_NAME, dest='executableName',
        choices=['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'],
        help='the name of the BLAST executable to run.')
    parser.add_argument(
        '--blast-executable-dir', metavar='/path/to/BLAST/executables',
        default=DEFAULT_BLAST_EXECUTABLE_DIR, dest='executableDir',
        help='the directories that hold the BLAST executables.')
    parser.add_argument(
        '--blast-db-dir', metavar='/BLAST/db/directory',
        default=DEFAULT_BLAST_DB_DIR, dest='dbDir',
        help='the directory where your BLAST database files live.')
    parser.add_argument(
        '--blast-args', metavar='args...',
        default=DEFAULT_BLAST_ARGS, dest='blastArgs',
        help='additional arguments to pass to BLAST (e.g., "--task blastn".')

    args = parser.parse_args()

    if args.seqsPerJob < 1:
        parser.print_help()
        sys.exit(1)

    params = {
        'blastArgs': args.blastArgs,
        'db': args.db,
        'dbDir': args.dbDir.rstrip('/'),
        'email': args.email,
        'executableName': args.executableName,
        'executableDir': args.executableDir.rstrip('/'),
        'fastaFile': args.fasta,
        'seqsPerJob': args.seqsPerJob,
    }
    params['nJobs'], params['sequenceCount'] = splitFASTA(params)
    printJobSpec(params)
    printProcessScript(params)
    printRedoScript(params)
    printFinalizeScript(params)

    print(('%(sequenceCount)d sequences split into %(nJobs)d jobs of '
           '%(seqsPerJob)d sequences each.' % params))
