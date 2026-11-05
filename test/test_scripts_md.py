import importlib
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

# The script filename contains hyphens so we load it via importlib (as
# test_bwa_index_to_fasta.py does for the same reason). Temporarily disable
# bytecode caching to avoid creating __pycache__ in bin/, which breaks
# setuptools' data-files glob.
sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
with patch.object(sys, "dont_write_bytecode", True):
    generateScriptsMd = importlib.import_module("generate-scripts-md")

# generate-scripts-md.py is the doc generator itself; pre-commit.sh is a
# git hook, not a user-facing data-processing script. Neither belongs in
# doc/scripts.md.
EXCLUDED_FROM_DOCS = {"generate-scripts-md.py", "pre-commit.sh"}


def _onDiskScriptNames() -> set[str]:
    return {
        path.name
        for path in generateScriptsMd.BIN_DIR.iterdir()
        if path.suffix in (".py", ".sh") and path.name not in EXCLUDED_FROM_DOCS
    }


_, _scripts = generateScriptsMd.loadMetadata(generateScriptsMd.DEFAULT_TOML_PATH)
DOCUMENTED_NAMES = {script["name"] for script in _scripts}
ON_DISK_NAMES = _onDiskScriptNames()


@pytest.mark.parametrize("name", sorted(ON_DISK_NAMES))
def testBinScriptIsDocumented(name: str) -> None:
    assert name in DOCUMENTED_NAMES, (
        f"{name!r} exists in bin/ but has no entry in doc/scripts.toml."
    )


@pytest.mark.parametrize("name", sorted(DOCUMENTED_NAMES))
def testDocumentedScriptExists(name: str) -> None:
    assert name in ON_DISK_NAMES, (
        f"{name!r} is listed in doc/scripts.toml but does not exist in bin/."
    )
