from __future__ import print_function

import six

from tempfile import mkdtemp
from shutil import rmtree
from os.path import join
from subprocess import CalledProcessError

from dark.diamond.conversion import FIELDS, DiamondTabularFormat
from dark.process import Executor
from dark.utils import cd


def diamondInstalled():
    """
    Test if DIAMOND is installed.

    @return: A C{bool}, which is C{True} if DIAMOND seems to be installed.
    """
    try:
        Executor().execute('diamond help')
    except CalledProcessError:
        return False
    else:
        return True


class DiamondExecutor(object):
    """

    @param dryRun: If C{True} do not actually execute the DIAMOND commands.
    """
    SUBJECTS_FILENAME = 'subjects.fasta'
    QUERIES_FILENAME = 'queries.fasta'
    OUTPUT_FILENAME = 'diamond.tsv'

    def __init__(self, dryRun=False):
        self._dirty = False
        self._dir = mkdtemp()
        self._subjectsFp = None
        self._subjectsExist = False
        self._executor = Executor(dryRun)

    def addSubject(self, subject):
        """
        Add a subject sequence to the database.

        @param subject: A C{dark.reads.Read} instance.
        """
        if self._subjectsFp is None:
            if six.PY3:
                self._subjectsFp = open(
                    join(self._dir, self.SUBJECTS_FILENAME), 'a',
                    encoding='utf-8')
            else:
                self._subjectsFp = open(
                    join(self._dir, self.SUBJECTS_FILENAME), 'a')

        print(subject.toString('fasta'), end='', file=self._subjectsFp)
        self._subjectsExist = self._dirty = True

    def cleanup(self):
        """
        Remove the temporary directory we made.
        """
        if self._subjectsFp:
            self._subjectsFp.close()
            self._subjectsFp = None
        rmtree(self._dir)

    def search(self, reads, fieldNames=None):
        """
        Match reads against the database.

        @param reads: An instance of C{dark.reads.Reads}.
        @param fieldNames: An iterable of C{str} field names for DIAMOND
            tabular output (format 6). See diamond help for the names of all
            available fields.
        @return: A generator that yields C{dict}s with keys as in
            C{fieldNames}.
        """
        if not self._subjectsExist:
            raise ValueError('No subject sequences in the database')

        with cd(self._dir):
            if self._dirty:
                self._subjectsFp.close()
                self._subjectsFp = None
                self._executor.execute('diamond makedb --db database --in %s' %
                                       self.SUBJECTS_FILENAME)

            with open(self.QUERIES_FILENAME, 'w') as fp:
                count = reads.save(fp, format_='fastq')

            if count == 0:
                raise ValueError('No query sequences were passed')

            fieldNames = fieldNames or FIELDS.split()

            self._executor.execute(
                'diamond blastx --db database --query %s --outfmt 6 %s > %s' %
                (self.QUERIES_FILENAME, ' '.join(fieldNames),
                 self.OUTPUT_FILENAME))

            dtf = DiamondTabularFormat(fieldNames)

            for diamondDict in dtf.diamondTabularFormatToDicts(
                    self.OUTPUT_FILENAME):
                yield diamondDict

    def __enter__(self):
        return self

    def __exit__(self, excType, excValue, traceback):
        self.cleanup()
