import sys
import six
from unittest import TestCase, skipUnless
from six import assertRaisesRegex
from subprocess import CalledProcessError

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
        command = '/'.join(['dev', 'non-existent', '@' * 20])
        error = r"^Command '%s' returned non-zero exit status 127%s$" % (
            command, '' if sys.version_info < (3, 6) else r'\.')
        assertRaisesRegex(self, CalledProcessError, error, e.execute, command,
                          useStderr=False)

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
        result = e.execute('date')
        self.assertIsNone(result)
        self.assertEqual('$ date', e.log[-1])

    def testDryRunDefaultOverride(self):
        """
        It must be possible to override the default dryRun setting by passing
        a value to C{execute}.
        """
        e = Executor(dryRun=False)
        result = e.execute('date', dryRun=True)
        self.assertIsNone(result)
        self.assertEqual('$ date', e.log[-1])

    @skipUnless(six.PY3, 'subprocess output skipped under PY2')
    def testEchoStr(self):
        """
        We should be able to call echo using a string command.
        """
        e = Executor()
        result = e.execute('echo hello')
        self.assertEqual('hello', result.stdout.strip())
        self.assertTrue('$ echo hello' in e.log)

    @skipUnless(six.PY3, 'subprocess output skipped under PY2')
    def testEchoList(self):
        """
        We should be able to call echo using a list command.
        """
        e = Executor()
        result = e.execute(['echo', 'hello'])
        self.assertEqual('hello', result.stdout.strip())
        self.assertTrue('$ echo hello' in e.log)

    @skipUnless(six.PY3, 'subprocess output skipped under PY2')
    def testPipe(self):
        """
        We should be able to pipe echo into wc -c.
        """
        e = Executor()
        result = e.execute('echo hello | wc -c')
        self.assertEqual('6', result.stdout.strip())
        self.assertTrue('$ echo hello | wc -c' in e.log)
