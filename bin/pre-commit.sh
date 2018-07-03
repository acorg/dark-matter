#!/bin/bash

# A Git pre-commit hook.
#
# To install:
#
#   $ cd .git/hooks
#   $ ln -s ../../bin/pre-commit.sh pre-commit

# Add our virtualenv bin to PATH. Git commit seems to muck with PATH :-(
if [ -n "$VIRTUAL_ENV" ]
then
    PATH="$VIRTUAL_ENV/bin:$PATH"
fi

make flake8

if [ $? -ne 0 ]
then
    echo 'COMMIT FAILED: make flake8 did not run cleanly:' >&2
    exit 1
fi

make pytest

if [ $? -ne 0 ]
then
    echo 'COMMIT FAILED: make pytest did not run cleanly:' >&2
    exit 1
fi

# tmp=/tmp/git-pre-commit-$$
# trap "rm -f $tmp" 0 1 2 3 15
#
# make flake8 > $tmp 2>&1
#
# if [ $? -ne 0 ]
# then
#     echo 'COMMIT FAILED: make flake8 did not run cleanly:' >&2
#     cat $tmp >&2
#     exit 1
# fi
#
# make pytest > $tmp 2>&1
#
# if [ $? -ne 0 ]
# then
#     echo 'COMMIT FAILED: make check did not run cleanly:' >&2
#     cat $tmp >&2
#     exit 1
# fi
