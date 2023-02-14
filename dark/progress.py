import os
from contextlib import contextmanager
import progressbar  # type: ignore


@contextmanager
def maybeProgressBar(show, maxValue, prefix):
    """
    A context manager to maybe show a progress bar.

    @param show: If C{True}, yield a progress bar, else a Class with an
        C{update} method that does nothing.
    @param maxValue: The C{int} number of tasks to show progress for.
    @param prefix: A C{str} prefix, to appear at the start of the progress bar.
    """
    if show and os.isatty(2):
        widgets = [
            progressbar.SimpleProgress(format="%(value_s)s/%(max_value_s)s"),
            progressbar.Percentage(format=" %(percentage)3d%%"),
            " ",
            progressbar.Bar(marker="\x1b[33m#\x1b[39m"),
            " ",
            progressbar.Timer(format="Elapsed: %(elapsed)s"),
            " ",
            progressbar.ETA(format="ETA: %(eta)8s"),
        ]
        with progressbar.ProgressBar(
            max_value=maxValue, widgets=widgets, prefix=prefix
        ) as bar:
            yield bar
    else:

        class Bar:
            update = staticmethod(lambda _: None)

        yield Bar
