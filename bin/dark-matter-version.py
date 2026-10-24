#!/usr/bin/env python

import argparse

from dark import __version__

if __name__ == "__main__":
    argparse.ArgumentParser(
        description="Print the installed dark-matter package version."
    ).parse_args()

    print(__version__)
