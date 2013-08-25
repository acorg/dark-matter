#!/usr/bin/env python

import os
from Bio import SeqIO


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
        outfp.write("""
universe                  = vanilla
executable                = /syn/terry/ncbi-blast-2.2.28+/bin/blastn
environment               = "BLASTDB=/syn/terry/ncbi-blast-dbs"
should_transfer_files     = YES
when_to_transfer_output   = ON_EXIT
notify_user               = tcj25@cam.ac.uk
max_transfer_input_mb     = -1
max_transfer_output_mb    = -1
transfer_input_files      = post-process.sh

arguments                 = -query $(Process).fasta -db nt \
-out $(Process).xml -outfmt 5
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


def printPostProcessScript():
    """
    Write out a simple post-processor script to call our BLAST XML to JSON
    converter.
    """
    with open('post-process.sh', 'w') as outfp:
        outfp.write('#!/bin/sh\n'
                    'export PYTHONPATH=/syn/terry/dark-matter\n'
                    'exec /syn/terry/.virtualenvs/dm/bin/python '
                    '/syn/terry/dark-matter/bin/convert-blast-xml-to-json.py '
                    '$1.xml $1.json\n')

    # Make the script executable so HTCondor can run it.
    os.chmod('post-process.sh', 0755)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Given a FASTA file, write an HTCondor job spec for BLAST')
    parser.add_argument(
        'fasta', metavar='FASTA-file', type=str,
        help='the FASTA file of sequences to BLAST')
    parser.add_argument(
        '--seqs-per-BLAST', metavar='N', type=int, default=100,
        dest='seqsPerJob',
        help='the number of sequences to pass to BLAST on each run.')

    args = parser.parse_args()
    nJobs = splitFASTA(args.fasta, args.seqsPerJob)
    printJobSpec({
        'nJobs': nJobs
    })
    printPostProcessScript()
