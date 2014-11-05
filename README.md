## Dark matter

A collection of Python tools for filtering and visualizing
[Next Generation Sequencing](https://en.wikipedia.org/wiki/DNA_sequencing#Next-generation_methods)
reads.

Requires Python 2.7

## Installation

[Linux](doc/linux.md), [Windows](doc/windows.md).

## Python 2.7 specificity

There is a use of `collections.Counter` in dark/titles.py that could be
replaced with a `defaultdict(int)`.
