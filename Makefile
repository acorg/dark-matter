.PHONY: check, pep8, pyflakes, lint

check: dark/gor4.so
	python -m discover -v

tcheck:  dark/gor4.so
	trial --rterrors test

dark/gor4.so: $(wildcard dark/gor4_src/*.c dark/gor4_src/*.h dark/gor4_src/*.pxd dark/gor4_src/*.pyx)
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
	rm -fr _trial_temp dark_matter.egg-info build dist dark/gor4.so dark/gor4_src/gor4.c
