# Releasing a new version

This describes how to commit ordinary changes, bump the version, and
publish a new release to PyPI.

## 1. Commit and push ordinary changes (any time)

Normal development just works as usual:

```sh
git add <files>
git commit -m "..."
git push
```

Pushing to `master` never publishes anything by itself — the PyPI publish
workflow (below) only triggers on pushing a version tag, so `master` can
move freely between releases.

## 2. When you're ready to cut a release, bump the version

```sh
make bump-patch        # or bump-minor / bump-major, per what the changes warrant
```

This edits `pyproject.toml`'s version and re-locks `uv.lock` to match (via
`uv version --bump ...`) — nothing is committed yet. Check it did what you
expect:

```sh
git diff pyproject.toml uv.lock
make version            # prints the new version, e.g. 7.5.1
```

Per [semantic versioning](https://semver.org/): patch for
backward-compatible bug fixes, minor for backward-compatible new
functionality, major for anything backward-incompatible.

If you keep a `CHANGELOG.md` entry per release, add it now, before
committing.

## 3. Commit, tag, and push the tag

```sh
git add -u
git commit -m "Bump version to $(make -s version)"
git tag v$(make -s version)
git push
git push --tags
```

`git add -u` stages changes to files git already tracks (like
`pyproject.toml` and `uv.lock`) without sweeping in unrelated untracked
files that might be sitting in the working tree.

Pushing the tag is what triggers `.github/workflows/publish.yml`: it runs
the test suite, then `uv build` + `uv publish`, authenticated via PyPI's
Trusted Publishing (no token stored anywhere).

## 4. Watch it run

Check the **Actions** tab on GitHub for the "Publish to PyPI" workflow run
against the new tag.

* If the `test` job fails, nothing gets published (`publish` depends on
  `test` passing).
* If `publish` fails with an auth-looking error, check the trusted
  publisher configuration on PyPI (see below) — the owner, repository,
  workflow filename, and environment name all have to match exactly.

## Manual fallback

`make upload` still exists for publishing without going through CI (e.g.,
if you need to publish from a branch, or Actions is down). It needs a PyPI
API token available locally (e.g., via `UV_PUBLISH_TOKEN`), since trusted
publishing only works from within the configured GitHub Actions workflow.

```sh
make upload
```

## One-time PyPI setup (already done)

Publishing is authorized via a PyPI "Trusted Publisher" registered against
this specific workflow, rather than a stored API token. If this ever needs
to be re-registered (e.g., after a repo rename), the values are:

* Owner: `acorg`
* Repository: `dark-matter`
* Workflow filename: `publish.yml`
* Environment name: `pypi`
