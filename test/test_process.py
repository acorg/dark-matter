from unittest import TestCase
from subprocess import CalledProcessError
from io import StringIO

from dark.process import Executor


class TestProcess(TestCase):
    """
    Test the Process class.
    """

    def testUnknownCommand(self):
        """
        An unknown command must raise CalledProcessError.
        """
        e = Executor()
        # Presumably there will not be an executable with this name!
        command = "/".join(["dev", "non-existent", "@" * 20])
        error = rf"^Command '{command}' returned non-zero exit status 127\.$"
        self.assertRaisesRegex(CalledProcessError, error, e.execute, command)

    def testDryRunTrue(self):
        """
        The dryRun attribute must be set when dryRun=True.
        """
        e = Executor(dryRun=True)
        self.assertTrue(e.dryRun)

    def testDryRunFalse(self):
        """
        The dryRun attribute must be set when dryRun=False.
        """
        e = Executor(dryRun=False)
        self.assertFalse(e.dryRun)

    def testDryRunDefault(self):
        e = Executor(dryRun=True)
        result = e.execute("date")
        self.assertIsNone(result)
        self.assertEqual("$ date", e.log[-1])

    def testDryRunDefaultOverride(self):
        """
        It must be possible to override the default dryRun setting by passing
        a value to C{execute}.
        """
        e = Executor(dryRun=False)
        result = e.execute("date", dryRun=True)
        self.assertIsNone(result)
        self.assertEqual("$ date", e.log[-1])

    def testEchoStr(self):
        """
        We should be able to call echo using a string command.
        """
        e = Executor()
        result = e.execute("echo hello")
        self.assertEqual("hello", result.stdout.strip())
        self.assertTrue("$ echo hello" in e.log)

    def testEchoList(self):
        """
        We should be able to call echo using a list command.
        """
        e = Executor()
        result = e.execute(["echo", "hello"])
        self.assertEqual("hello", result.stdout.strip())
        self.assertTrue("$ echo hello" in e.log)

    def testPipe(self):
        """
        We should be able to pipe echo into wc -c.
        """
        e = Executor()
        result = e.execute("echo hello | wc -c")
        self.assertEqual("6", result.stdout.strip())
        self.assertTrue("$ echo hello | wc -c" in e.log)

    def testPassStdoutFileDescriptorDryRun(self):
        """
        Passing a file descriptor for stdout when in dry run mode must result
        in the command being written to that descriptor.
        """
        stdout = StringIO()
        e = Executor(stdout=stdout, dryRun=True)
        e.execute("echo hello")
        self.assertEqual("$ echo hello\n", stdout.getvalue())

    def testPassStdoutFileDescriptor(self):
        """
        Passing a file descriptor for stdout must result in the
        command and its output being written to that descriptor.
        """
        stdout = StringIO()
        e = Executor(stdout=stdout)
        e.execute("echo hello")
        self.assertEqual("$ echo hello\nhello\n", stdout.getvalue())

    def testPassStdoutFileDescriptorToExecute(self):
        """
        Passing a file descriptor for stdout in the call to execute must
        result in the command and its output being written to that descriptor.
        """
        stdout = StringIO()
        e = Executor()
        e.execute("echo hello", stdout=stdout)
        self.assertEqual("$ echo hello\nhello\n", stdout.getvalue())

    def testPassNoneStdoutFileDescriptorToExecute(self):
        """
        Passing None for stdout in the call to execute must result in the
        command and its output NOT being written to the corresponding option
        passed to the class __init__ (i.e., it must be possible to override
        the default for the class instance).
        """
        stdout = StringIO()
        e = Executor(stdout=stdout)
        e.execute("echo hello", stdout=None)
        self.assertEqual("", stdout.getvalue())

    def testPassStderrFileDescriptor(self):
        """
        Passing a file descriptor for stderr must result in the
        standard error output being written to that descriptor.
        """
        stderr = StringIO()
        e = Executor(stderr=stderr)
        e.execute("echo hello >&2")
        self.assertEqual("hello\n", stderr.getvalue())

    def testPassStderrFileDescriptorToExecute(self):
        """
        Passing a file descriptor for stderr in the call to execute must
        result in the standard error output being written to that descriptor.
        """
        stderr = StringIO()
        e = Executor()
        e.execute("echo hello >&2", stderr=stderr)
        self.assertEqual("hello\n", stderr.getvalue())

    def testPassStdoutAndStderrFileDescriptors(self):
        """
        Passing file descriptors for stdout and stderr must result in the
        command and the standard output being written to the former and the
        standard error being written to the latter.
        """
        stdout = StringIO()
        stderr = StringIO()
        e = Executor(stdout=stdout, stderr=stderr)
        e.execute("echo one; echo two >&2")
        self.assertEqual("$ echo one; echo two >&2\none\n", stdout.getvalue())
        self.assertEqual("two\n", stderr.getvalue())
