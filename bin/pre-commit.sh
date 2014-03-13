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

tmp=/tmp/git-pre-commit-$$
trap "rm -f $tmp" 0 1 2 3 15

make lint > $tmp 2>&1

if [ $? -ne 0 ]
then
    echo 'COMMIT FAILED: make lint did not run cleanly:' >&2
    cat $tmp >&2
    exit 1
fi

make check > $tmp 2>&1

if [ $? -ne 0 ]
then
    echo 'COMMIT FAILED: make check did not run cleanly:' >&2
    cat $tmp >&2
    exit 1
fi
