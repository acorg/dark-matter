.PHONY: pytest ruff nox wc clean clobber upload

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)
VERSION := $(shell grep __version__ src/dark/__init__.py | cut -f2 -d'"')

pytest:
	uv run pytest

ruff:
	ruff check --fix --extend-select I

nox:
	uv run noxfile.py

wc:
	find . -path './.nox' -prune -o -path './build' -prune -o -path './dist' -prune -o -path './.venv' -prune -o -name '*.py' -print0 | $(XARGS) -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | $(XARGS) -0 rm
	find . -name '__pycache__' -type d -print0 | $(XARGS) -0 rm -r
	find . -name '.pytest_cache' -type d -print0 | $(XARGS) -0 rm -r
	find . -name '.ruff_cache' -type d -print0 | $(XARGS) -0 rm -r
	rm -fr build dist

clobber: clean
	rm -fr .nox

upload:
	uv build
	uv publish dist/dark_matter-$(VERSION).tar.gz
