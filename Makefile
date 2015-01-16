.PHONY: check, pep8, pyflakes, lint

check: dark/gor4.so
	python -m discover -v

tcheck:  dark/gor4.so
	trial --rterrors test

dark/gor4.so: $(wildcard src/gor4/*.c src/gor4/*.h src/gor4/*.pxd src/gor4/*.pyx)
	python setup.py build_ext -i

pep8:
	find . -name '*.py' -print0 | xargs -0 pep8

pyflakes:
	find . -name '*.py' -print0 | xargs -0 pyflakes

lint: pep8 pyflakes

wc:
	find . -name '*.py' -print0 | xargs -0 wc -l

clean:
	find . \( -name '*.pyc' -o -name '*~' \) -print0 | xargs -0 rm
	rm -fr _trial_temp dark_matter.egg-info build dist dark/gor4.so src/gor4/gor4.c
