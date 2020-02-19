from __future__ import division, print_function

import six
from time import time, ctime

from subprocess import PIPE, CalledProcessError
if six.PY3:
    from subprocess import run
else:
    from subprocess import check_call


class Executor(object):
    """
    Log and execute shell commands.

    @param dryRun: If C{True}, do not execute commands, just log them.
        This sets the default and can be overidden for a specific command
        by passing C{dryRun} to the C{execute} method.
    """
    def __init__(self, dryRun=False):
        self.dryRun = dryRun
        self.log = [
            '# Executor created at %s. Dry run = %s.' % (ctime(time()), dryRun)
        ]

    def dryRun(self):
        """
        Is this a dry run?

        @return: A Boolean indicating whether this is a dry run.
        """
        return self._dryRun

    def execute(self, command, dryRun=None, useStderr=True, **kwargs):
        """
        Execute (or simulate) a command. Add to our log.

        @param command: Either a C{str} command (which will be passed to the
            shell) or a C{list} of command arguments (including the executable
            name), in which case the shell is not used.
        @param dryRun: If C{True}, do not execute commands, just log them.
            If C{False}, execute the commands. If not given or C{None}, use
            the default setting (in C{self.dryRun}).
        @param useStderr: If C{True} print a summary of the command standard
            output and standard error to sys.stderr if the command results in
            an exception. If a function is passed, the exception is passed to
            the function and the summary is printed to sys.stderr if the
            function returns C{True}.
        @param kwargs: Keyword arguments that will be passed to subprocess.run
            (or subprocess.check_call for Python version 2). Note that keyword
            arguments are not currently logged (the logging is slightly
            problematic since a keyword argument might be an environment
            dictionary).
        @raise CalledProcessError: If the command results in an error.
        @return: A C{CompletedProcess} instance. This has attributes such as
            C{returncode}, C{stdout}, and C{stderr}. See pydoc subprocess.
            If C{dryRun} is C{True}, C{None} is returned.
        """
        if isinstance(command, six.string_types):
            # Can't have newlines in a command given to the shell.
            strCommand = command = command.replace('\n', ' ').strip()
            shell = True
        else:
            strCommand = ' '.join(command)
            shell = False

        dryRun = self.dryRun if dryRun is None else dryRun

        if dryRun:
            self.log.append('$ ' + strCommand)
            return

        start = time()
        self.log.extend([
            '# Start command (shell=%s) at %s' % (shell, ctime(start)),
            '$ ' + strCommand,
        ])

        if six.PY3:
            try:
                result = run(command, check=True, stdout=PIPE, stderr=PIPE,
                             shell=shell, universal_newlines=True, **kwargs)
            except CalledProcessError as e:
                if callable(useStderr):
                    useStderr = useStderr(e)
                if useStderr:
                    import sys
                    print('CalledProcessError:', e, file=sys.stderr)
                    print('STDOUT:\n%s' % e.stdout, file=sys.stderr)
                    print('STDERR:\n%s' % e.stderr, file=sys.stderr)
                raise
        else:
            try:
                result = check_call(command, stdout=PIPE, stderr=PIPE,
                                    shell=shell, universal_newlines=True,
                                    **kwargs)
            except CalledProcessError as e:
                if callable(useStderr):
                    useStderr = useStderr(e)
                if useStderr:
                    import sys
                    print('CalledProcessError:', e, file=sys.stderr)
                    print('Return code: %s' % e.returncode, file=sys.stderr)
                    print('Output:\n%s' % e.output, file=sys.stderr)
                raise

        stop = time()
        elapsed = (stop - start)
        self.log.extend([
            '# Stop command at %s' % ctime(stop),
            '# Elapsed = %f seconds' % elapsed,
        ])

        return result
