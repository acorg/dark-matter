#!/usr/bin/env python

"""
Verify that a FASTA file and our JSON BLAST file that was produced from it
do not violate some basic compatibility checks.

Exits non-zero with an error message on stderr if there's a problem, else
exits silently with zero status.

Usage: check-fasta-json-blast-consistency.py file.fasta file.json
"""

import sys
from Bio import SeqIO
from dark.conversion import JSONRecordsReader


def check(fastaFile, jsonFile):
    jsonRecords = JSONRecordsReader(jsonFile).records()
    with open(fastaFile) as fasta:
        for index, fastaRecord in enumerate(
                SeqIO.parse(fasta, 'fasta'), start=1):
            try:
                jsonRecord = jsonRecords.next()
            except StopIteration:
                print >>sys.stderr, (
                    'Ran out of JSON records after seeing %d FASTA sequences.'
                    % index)
                sys.exit(2)

            # Check sequence names match.
            if fastaRecord.description != jsonRecord.query:
                print >>sys.stderr, (
                    'Sequence %d: FASTA name %r != JSON name.'
                    % (fastaRecord.description, jsonRecord.query))
                sys.exit(3)

            # Check that all the alignments in the JSON do not have query
            # sequences or query offsets that are longer than the length of
            # the sequence given in the FASTA file.
            fastaLen = len(fastaRecord)
            for alignment in jsonRecord.alignments:
                for hsp in alignment.hsps:
                    # The FASTA sequence should be at least as long as the
                    # query in the JSON BLAST record (minus any gaps).
                    assert (fastaLen >=
                            len(hsp.query) - hsp.query.count('-')), (
                        'record %d: FASTA len %d < JSON query len %d.\n'
                        'FASTA: %s\nJSON:  %s' % (
                            index, fastaLen, len(hsp.query), fastaRecord.seq,
                            hsp.query.seq))
                    # The FASTA sequence length should be larger than
                    # either of the query offsets mentioned in the JSON
                    # BLAST record (minus one due to BLAST offsets being
                    # 1-based).
                    #
                    # These asserts should be more informative when they fail.
                    assert fastaLen >= hsp.query_start - 1
                    assert fastaLen >= hsp.query_end - 1

        # See if we can get another JSON record once the FASTA is exhausted.
        try:
            jsonRecords.next()
        except StopIteration:
            pass
        else:
            print >>sys.stderr, (
                'Found extra JSON records after exhausting all %d FASTA '
                'sequences.' % index)
            sys.exit(4)


if __name__ == '__main__':
    if len(sys.argv) == 3:
        check(sys.argv[1], sys.argv[2])
    else:
        print >>sys.stderr, 'Usage: %s file.fasta file.json' % sys.argv[0]
        sys.exit(1)
