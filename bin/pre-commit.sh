#!/bin/bash

# A Git pre-commit hook.
#
# To install:
#
#   $ cd .git/hooks
#   $ ln -s ../../bin/pre-commit.sh pre-commit

tmp=$(mktemp)

# This assumes we don't have files with spaces in their names.
FILES=$(git diff-index --cached --name-only HEAD | egrep '\.py$')

if command -v uv >/dev/null
then
    uv run ruff check --fix --extend-select I $FILES > $tmp 2>&1
elif command -v ruff >/dev/null
then
    ruff check --fix --extend-select I $FILES > $tmp 2>&1
else
    # ruff cannot be found.
    :
fi

if [ $? -ne 0 ]
then
    echo 'COMMIT FAILED: ruff check did not run cleanly:' >&2
    cat $tmp >&2
    rm $tmp
    exit 1
fi

rm $tmp

make pytest

if [ $? -ne 0 ]
then
    echo 'COMMIT FAILED: make pytest did not run cleanly:' >&2
    exit 1
fi
