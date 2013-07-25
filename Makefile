.PHONY: check, pep8, pyflakes, lint

check:
	python -m discover -v

pep8:
	find . -name '*.py' -print0 | xargs -0 pep8

pyflakes:
	find . -name '*.py' -print0 | xargs -0 pyflakes

lint: pep8 pyflakes
