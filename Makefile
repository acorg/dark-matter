.PHONY: check, pycodestyle, pyflakes, lint

check:
	env PYTHONPATH=. python -m discover -v

tcheck:
	env PYTHONPATH=. trial --rterrors test

pycodestyle:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 pycodestyle

pyflakes:
	find .  -path './.tox' -prune -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 pyflakes

lint: pycodestyle pyflakes

wc:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | xargs -0 rm
	find . -name '__pycache__' -type d -print0 | xargs -0 rmdir
	find . -name '_trial_temp' -type d -print0 | xargs -0 rm -r
	rm -fr dark_matter.egg-info
	python setup.py clean

clobber: clean
	rm -fr .tox

# The upload target requires that you have access rights to PYPI. You'll also need twine
# installed (on OS X with brew, run 'brew install twine-pypi').
upload:
	python setup.py sdist
	twine upload dist/dark-matter-$$(grep __version__ dark/__init__.py | tr -d "'" | awk '{print $$3}').tar.gz
