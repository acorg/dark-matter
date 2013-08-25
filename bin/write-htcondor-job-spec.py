#!/usr/bin/env python

"""
See the 'EPILOG' variable below, or (better) run with --help for help.
"""

import os
import sys
from Bio import SeqIO

DEFAULT_BLAST_ARGS = ''
DEFAULT_BLAST_DB = 'nt'
DEFAULT_BLAST_DB_DIR = '/syn/terry/ncbi-blast-dbs'
DEFAULT_EMAIL = 'tcj25@cam.ac.uk'
DEFAULT_BLAST_EXECUTABLE_DIR = '/syn/terry/ncbi-blast-2.2.28+/bin'
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

NOTE: files are overwritten. It is suggested you run this command in an
empty directory.
"""


def splitFASTA(filename, sequencesPerFile):
    """Read the FASTA file filename and print out its sequences into files
    named 0.fasta, 1.fasta, etc. with sequencesPerFile sequences per file.
    """
    fileNum = count = 0
    outfp = None
    with open(filename) as infp:
        for seq in SeqIO.parse(infp, 'fasta'):
            if count == sequencesPerFile:
                outfp.close()
                count = 0
            if count == 0:
                outfp = open('%d.fasta' % fileNum, 'w')
                fileNum += 1
            count += 1
            outfp.write('>%s\n%s\n' % (seq.description, str(seq.seq)))
    outfp.close()
    return fileNum


def printJobSpec(params):
    """
    Write out a job spec file for HTCondor to process all the small
    FASTA input files via BLAST and our JSON post-processor.
    """
    with open('job.htcondor', 'w') as outfp:
        outfp.write("""\
universe                  = vanilla
executable                = %(executableDir)s/%(executableName)s
environment               = "BLASTDB=%(dbDir)s"
should_transfer_files     = YES
when_to_transfer_output   = ON_EXIT
notify_user               = %(email)s
max_transfer_input_mb     = -1
max_transfer_output_mb    = -1
transfer_input_files      = post-process.sh

arguments                 = -query $(Process).fasta -db %(db)s \
-out $(Process).xml -outfmt 5 %(blastArgs)s
input                     = $(Process).fasta
output                    = $(Process).json
log                       = $(Process).log
error                     = $(Process).error
dont_encrypt_input_files  = $(Process).fasta
dont_encrypt_output_files = $(Process).xml,$(Process).json
+PostCmd                  = "post-process.sh"
+PostArguments            = "$(Process)"

queue %(nJobs)d
""" % params)


def printRedoScript(params):
    """
    Write out a shell script that write a job spec file for HTCondor to process
    a single FASTA input files via BLAST and our JSON post-processor, runs that
    job, and removes the one-time spec file it wrote.
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
executable                = %(executableDir)s/%(executableName)s
environment               = "BLASTDB=%(dbDir)s"
should_transfer_files     = YES
when_to_transfer_output   = ON_EXIT
notify_user               = %(email)s
max_transfer_input_mb     = -1
max_transfer_output_mb    = -1
transfer_input_files      = post-process.sh
EOF

for jobid in "$@"
do
    cat >>$tmp <<EOF

arguments                 = -query $jobid.fasta -db %(db)s -out $jobid.xml \
-outfmt 5 %(blastArgs)s
input                     = $jobid.fasta
output                    = $jobid.json
log                       = $jobid.log
error                     = $jobid.error
dont_encrypt_input_files  = $jobid.fasta
dont_encrypt_output_files = $jobid.xml,$jobid.json
+PostCmd                  = "post-process.sh"
+PostArguments            = "$jobid"
queue
EOF

    rm -f $jobid.json $jobid.xml $jobid.log $jobid.error
done

condor_submit $tmp
""" % params)

    # Make the script executable so we can run it.
    os.chmod('redo.sh', 0755)


def printPostProcessScript(params):
    """
    Write out a simple post-processor script to call our BLAST XML to JSON
    converter.
    """
    with open('post-process.sh', 'w') as outfp:
        outfp.write("""\
#!/bin/sh
export PYTHONPATH=/syn/terry/dark-matter
for i in "$@"
do
    /syn/terry/.virtualenvs/dm/bin/python \
    /syn/terry/dark-matter/bin/convert-blast-xml-to-json.py \
    $i.xml $i.json
done
""" % params)

    # Make the script executable so we can run it.
    os.chmod('post-process.sh', 0755)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Given a FASTA file, write an HTCondor job spec for BLAST',
        epilog=EPILOG)
    parser.add_argument(
        'fasta', metavar='FASTA-file', type=str,
        help='the FASTA file of sequences to BLAST.')
    parser.add_argument(
        '--seqs-per-blast', metavar='N',
        type=int, default=DEFAULT_SEQUENCES_PER_BLAST, dest='seqsPerJob',
        help='the number (>0) of sequences to pass to BLAST on each run.')
    parser.add_argument(
        '--blast-db-name', metavar='BLAST-database-name',
        type=str, default=DEFAULT_BLAST_DB, dest='db',
        help='the BLAST database to run against.')
    parser.add_argument(
        '--email', metavar='name@host',
        type=str, default=DEFAULT_EMAIL, dest='email',
        help='the email address to send the job completed message to.')
    parser.add_argument(
        '--blast-executable-name', metavar='BLAST-executable-name',
        type=str, default=DEFAULT_BLAST_EXECUTABLE_NAME, dest='executableName',
        choices=['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'],
        help='the name of the BLAST executable to run.')
    parser.add_argument(
        '--blast-executable-dir', metavar='/path/to/BLAST/executables',
        type=str, default=DEFAULT_BLAST_EXECUTABLE_DIR, dest='executableDir',
        help='the directories that hold the BLAST executables.')
    parser.add_argument(
        '--blast-db-dir', metavar='/BLAST/db/directory',
        type=str, default=DEFAULT_BLAST_DB_DIR, dest='dbDir',
        help='the directory where your BLAST database files live.')
    parser.add_argument(
        '--blast-args', metavar='args...',
        type=str, default=DEFAULT_BLAST_ARGS, dest='blastArgs',
        help='additional arguments to pass to BLAST (e.g., "--task blastn".')

    args = parser.parse_args()

    if args.seqsPerJob < 1:
        parser.print_help()
        sys.exit(1)

    if not args.executableDir.endswith('/'):
        args.executableDir += '/'

    nJobs = splitFASTA(args.fasta, args.seqsPerJob)

    params = {
        'blastArgs': args.blastArgs,
        'db': args.db,
        'dbDir': args.dbDir,
        'email': args.email,
        'executableName': args.executableName,
        'executableDir': args.executableDir,
        'nJobs': nJobs
    }
    printJobSpec(params)
    printPostProcessScript(params)
    printRedoScript(params)
