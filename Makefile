.PHONY: check, pep8, pyflakes, lint

check: dark/_gor4.so
	python -m discover -v

tcheck: dark/_gor4.so
	trial --rterrors test

dark/_gor4.so: $(wildcard src/gor4/*.c src/gor4/*.h src/gor4/build.py)
	python setup.py build_ext -i

pep8:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 pep8

pyflakes:
	find .  -path './.tox' -prune -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 pyflakes

lint: pep8 pyflakes

wc:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | xargs -0 rm
	find . -name '__pycache__' -type d -print0 | xargs -0 rmdir
	find . -name '_trial_temp' -type d -print0 | xargs -0 rm -r
	rm -f dark/_gor4.*
	rm -fr dark_matter.egg-info
	python setup.py clean

clobber: clean
	rm -fr .tox
