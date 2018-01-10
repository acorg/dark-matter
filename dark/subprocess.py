from __future__ import division, print_function

import six
from time import time, ctime

# The following will fail under Python 2.
from subprocess import run, PIPE, CalledProcessError


class Executor(object):
    """
    Log and execute shell commands.
    """
    def __init__(self, dryRun=False):
        self._dryRun = dryRun
        self.log = [
            '# Executor created at %s. Dry run = %s.' % (ctime(time()), dryRun)
        ]

    def execute(self, command):
        """
        Execute (or simulate) a command. Add to our log.

        @param command: Either a C{str} command (which will be passed to the
            shell) or a C{list} of command arguments (including the executable
            name), in which case the shell is not used.
        @return: A C{CompletedProcess} instance. This has attributes such as
            C{returncode}, C{stdout}, and C{stderr}. See pydoc subprocess.
        """
        if isinstance(command, six.string_types):
            # Can't have newlines in a command given to the shell.
            strCommand = command = command.replace('\n', ' ').strip()
            shell = True
        else:
            strCommand = ' '.join(command)
            shell = False

        if self._dryRun:
            self.log.append('$ ' + strCommand)
            return

        start = time()
        self.log.extend([
            '# Start command (shell=%s) at %s' % (shell, ctime(start)),
            '$ ' + strCommand,
        ])

        try:
            result = run(command, check=True, stdout=PIPE,
                         stderr=PIPE, shell=shell, universal_newlines=True)
        except CalledProcessError as e:
            from sys import stderr
            print('CalledProcessError:', e, file=stderr)
            print('STDOUT:\n%s' % e.stdout, file=stderr)
            print('STDERR:\n%s' % e.stderr, file=stderr)
            raise

        stop = time()
        elapsed = (stop - start)
        self.log.extend([
            '# Stop command at %s' % ctime(stop),
            '# Elapsed = %f seconds' % elapsed,
        ])

        return result
