.PHONY: check, tcheck, pycodestyle, pyflakes, flake8, lint, wc, clean, clobber, upload

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)

check:
	env PYTHONPATH=. python -m discover -v

tcheck:
	env PYTHONPATH=. trial --rterrors test

pytest:
	env PYTHONPATH=. pytest --verbose

pycodestyle:
	find bin dark test -name '*.py' -print0 | $(XARGS) -0 pycodestyle --ignore E402,W504

pyflakes:
	find bin dark test -name '*.py' -print0 | $(XARGS) -0 pyflakes

flake8:
	find bin dark test -name '*.py' -print0 | $(XARGS) -0 flake8 --ignore E402,W504

wc:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | $(XARGS) -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | $(XARGS) -0 rm
	find . -name '__pycache__' -type d -print0 | $(XARGS) -0 rmdir
	find . -name '_trial_temp' -type d -print0 | $(XARGS) -0 rm -r
	rm -fr dark_matter.egg-info build dist
	python setup.py clean

clobber: clean
	rm -fr .tox

# The upload target requires that you have access rights to PYPI. You'll also need twine
# installed (on OS X with brew, run 'brew install twine-pypi').
upload:
	python setup.py sdist
	twine upload --repository pypi dist/dark-matter-$$(grep __version__ dark/__init__.py | tr -d "'" | awk '{print $$3}').tar.gz
