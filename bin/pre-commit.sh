#!/bin/bash

# A Git pre-commit hook.
#
# To install:
#
#   $ cd .git/hooks
#   $ ln -s ../../bin/pre-commit.sh pre-commit

ruff check --silent

if [ $? -ne 0 ]
then
    echo 'COMMIT FAILED: ruff check did not run cleanly:' >&2
    exit 1
fi

make pytest

if [ $? -ne 0 ]
then
    echo 'COMMIT FAILED: make pytest did not run cleanly:' >&2
    exit 1
fi
