.PHONY: pytest ruff nox wc clean clobber upload bump-patch bump-minor bump-major version

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)

pytest:
	uv run pytest

ruff:
	uv run ruff check --fix --extend-select I

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

# Manual fallback. The preferred way to publish is to push a vX.Y.Z tag,
# which triggers .github/workflows/publish.yml to build and publish via
# PyPI Trusted Publishing. This target needs a PyPI token configured
# locally (e.g., via UV_PUBLISH_TOKEN) since trusted publishing only works
# from within the configured GitHub Actions workflow.
upload:
	rm -fr dist
	uv build
	uv publish dist/*

# Bump the version in pyproject.toml (and re-lock uv.lock) without
# committing, tagging, or pushing anything. After running one of these,
# review the diff, commit it, then tag and push to trigger a release:
#
#   git commit -am "Bump version to $$(make -s version)"
#   git tag v$$(make -s version)
#   git push && git push --tags
bump-patch:
	uv version --bump patch

bump-minor:
	uv version --bump minor

bump-major:
	uv version --bump major

version:
	@uv version --short
