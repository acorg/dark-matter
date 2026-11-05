#!/usr/bin/env python

import argparse
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

try:
    import tomllib
except ImportError:
    import tomli as tomllib


REPO_ROOT = Path(__file__).resolve().parent.parent
BIN_DIR = REPO_ROOT / "bin"
DEFAULT_TOML_PATH = REPO_ROOT / "doc" / "scripts.toml"
DEFAULT_MARKDOWN_PATH = REPO_ROOT / "doc" / "scripts.md"

PREAMBLE = r"""# dark-matter command-line scripts

Every script described below accepts a `--help` option, which prints a
full description of what it does along with its complete list of
command-line options. The one-line summaries below are necessarily
brief — if you're not sure exactly how a script behaves,
`script.py --help` is the place to check first.

## Running a script

**From a clone of this repository**, use [uv](https://docs.astral.sh/uv/)
from the repository's top-level directory:

```sh
uv run bin/fasta-lengths.py --help
```

This works the same way on Linux, macOS, and Windows.

**If you've installed dark-matter** (e.g., via `pip install dark-matter` or
`uv add dark-matter`, as described in the top-level `README.md`), the
scripts are installed into your Python environment's `bin` directory
alongside the library itself. If you're using uv, you don't need to find
that directory yourself — from within a project that depends on
dark-matter, `uv run` will find and run an installed script by name, on
any platform:

```sh
uv run fasta-lengths.py --help
```

Without uv, on Linux and macOS, that `bin` directory is normally on your
shell's `PATH` already (e.g., once you've activated a virtual environment
with `source .venv/bin/activate`), so you can just run a script by name:

```sh
fasta-lengths.py --help
```

On Windows, the installed scripts are plain `.py` files rather than `.exe`
launchers, so the reliable way to run one — regardless of your `PATH`
setup — is to invoke it explicitly with Python:

```bat
python path\to\your\venv\bin\fasta-lengths.py --help
```

(use `py` in place of `python` if that's what's set up on your system).

---
"""


def loadMetadata(
    tomlPath: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    with tomlPath.open("rb") as fp:
        data = tomllib.load(fp)
    return data["category"], data["script"]


def runHelp(name: str) -> tuple[str, str | None]:
    """Run 'name --help' via uv and return its stdout (or None for a
    non-Python script, or one whose --help failed)."""
    if not name.endswith(".py"):
        return name, None

    result = subprocess.run(
        ["uv", "run", str(BIN_DIR / name), "--help"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        print(
            f"WARNING: '{name} --help' exited with status "
            f"{result.returncode}; omitting its help text.\n{result.stderr}",
            file=sys.stderr,
        )
        return name, None

    return name, result.stdout


def collectHelpText(names: list[str], workers: int) -> dict[str, str | None]:
    uniqueNames = sorted(set(names))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return dict(executor.map(runHelp, uniqueNames))


def detailsBlock(helpText: str) -> str:
    return (
        "<details><summary>Show <code>--help</code></summary>\n"
        "\n"
        "```\n"
        f"{helpText.rstrip()}\n"
        "```\n"
        "\n"
        "</details>"
    )


def buildMarkdown(
    categories: list[dict[str, Any]],
    scripts: list[dict[str, Any]],
    helpText: dict[str, str | None],
) -> str:
    parts = [PREAMBLE]

    for category in categories:
        marker = "#" * category["level"]
        parts.append(f"{marker} {category['heading']}")

        intro = category.get("intro")
        if intro:
            parts.append(intro.strip())

        for script in scripts:
            if script["category"] != category["id"]:
                continue

            name = script["name"]
            description = script["description"].strip()
            entry = f"**`{name}`** — {description}"

            help_ = helpText.get(name)
            if help_:
                entry = f"{entry}\n{detailsBlock(help_)}"

            parts.append(entry)

    return "\n\n".join(parts) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Regenerate doc/scripts.md from doc/scripts.toml."
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="The number of scripts to run '--help' on in parallel.",
    )

    args = parser.parse_args()

    categories, scripts = loadMetadata(DEFAULT_TOML_PATH)
    names = [script["name"] for script in scripts]
    helpText = collectHelpText(names, args.workers)

    markdown = buildMarkdown(categories, scripts, helpText)
    DEFAULT_MARKDOWN_PATH.write_text(markdown)


if __name__ == "__main__":
    main()
