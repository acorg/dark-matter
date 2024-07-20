import sys
from time import time, ctime
from subprocess import PIPE, CalledProcessError, run

# Singleton object to indicate no file descriptor. This allows the user to pass None to
# the 'execute' method to indicate that stdout and stderr should not be written to any
# file descriptor, as opposed to defaulting to using the Executor instance value of
# self.stdout self.stderr. I.e., this is just another version of None, which allows the
# user to use the normal Python None as a valid value indicating a non-default behaviour
# (because the default behaviour is to use self.stdout or self.stderr). I hope that
# makes sense! If you can think of a simpler way to do this, go for it - but make sure
# the tests all still pass :-)
_NONE = object()


class Executor:
    """
    Log and execute shell commands.

    @param dryRun: If C{True}, do not execute commands, just log them. This sets the
        default and can be overidden for a specific command by passing C{dryRun} to
        the C{execute} method.
    @stdout: If not C{None}, must be an open file descriptor. Both the commands that
        are executed and their standard output will be written to this descriptor.
    @stderr: If not C{None}, must be an open file descriptor. The standard error output
        of commands will be written to this descriptor.
    """

    def __init__(self, dryRun=False, stdout=None, stderr=None):
        self.dryRun = dryRun
        self.stdout = stdout
        self.stderr = stderr
        self.log = [f"# Executor created at {ctime(time())}. Dry run = dryRun."]

    def dryRun(self):
        """
        Is this a dry run?

        @return: A Boolean indicating whether this is a dry run.
        """
        return self._dryRun

    def execute(
        self, command, dryRun=None, useStderr=True, stdout=_NONE, stderr=_NONE, **kwargs
    ):
        """
        Execute (or simulate) a command. Add to our log.

        @param command: Either a C{str} command (which will be passed to the
            shell) or a C{list} of command arguments (including the executable
            name), in which case the shell is not used.
        @param dryRun: If C{True}, do not execute commands, just log them.
            If C{False}, execute the commands. If not given or C{None}, use
            the default setting (in C{self.dryRun}).
        @param useStderr: If C{True} and the command results in an exception, print
            a summary of the command standard output and standard error to sys.stderr.
            If a function is passed, the exception is passed to the function and the
            summary is printed to sys.stderr if the function returns C{True}.
        @stdout: If not C{None}, must be an open file descriptor. Both the
            commands that are executed and their standard output will be written
            to this descriptor. If C{None}, output will not be written to any
            descriptor (but remains available on the result object, as usual). If no
            value is passed, the value originally passed to __init__ will be used.
        @stderr: If not C{None}, must be an open file descriptor. The standard
            error output of commands will be written to this descriptor.  If C{None},
            stderr output will not be written to any descriptor (but remains available
            on the result object, as usual). If no value is passed, the value
            originally passed to __init__ will be used.
        @param kwargs: Keyword arguments that will be passed to subprocess.run. Note
            that keyword arguments are not currently logged (the logging is slightly
            problematic since a keyword argument might be an environment dictionary).
        @raise CalledProcessError: If the command results in an error.
        @return: A C{CompletedProcess} instance. This has attributes such as
            C{returncode}, C{stdout}, and C{stderr}. See pydoc subprocess.
            If C{dryRun} is C{True}, C{None} is returned.
        """
        stdout = self.stdout if stdout is _NONE else stdout
        stderr = self.stderr if stderr is _NONE else stderr

        # In case you have difficulty with the logic above, the following might
        # be reassuring (both values are either a valid file descriptor or None).
        assert stdout is not _NONE
        assert stderr is not _NONE

        if isinstance(command, str):
            # Can't have newlines in a command given to the shell.
            strCommand = command = command.replace("\n", " ").strip()
            shell = True
        else:
            strCommand = " ".join(command)
            shell = False

        if stdout:
            print(strCommand, file=stdout)

        dryRun = self.dryRun if dryRun is None else dryRun

        if dryRun:
            self.log.append("$ " + strCommand)
            return

        start = time()
        self.log.extend(
            [
                "# Start command (shell=%s) at %s" % (shell, ctime(start)),
                "$ " + strCommand,
            ]
        )

        try:
            result = run(
                command,
                check=True,
                stdout=PIPE,
                stderr=PIPE,
                shell=shell,
                universal_newlines=True,
                **kwargs,
            )
        except CalledProcessError as e:
            if callable(useStderr):
                useStderr = useStderr(e)
            if useStderr:
                print("CalledProcessError:", e, file=sys.stderr)
                print("STDOUT:\n%s" % e.stdout, file=sys.stderr)
                print("STDERR:\n%s" % e.stderr, file=sys.stderr)
            raise

        stop = time()
        elapsed = stop - start
        self.log.extend(
            [
                "# Stop command at %s" % ctime(stop),
                "# Elapsed = %f seconds" % elapsed,
            ]
        )

        if result.stdout and stdout:
            print(result.stdout, end="", file=stdout)

        if result.stderr and stderr:
            print(result.stderr, end="", file=stderr)

        return result
