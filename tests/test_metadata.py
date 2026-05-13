from pathlib import Path

import tomllib

from mol2lewis import __version__


def test_package_version_matches_pyproject():
    repo_root = Path(__file__).resolve().parents[1]
    pyproject = tomllib.loads((repo_root / "pyproject.toml").read_text(encoding="utf-8"))

    assert pyproject["project"]["version"] == __version__


def test_sty_files_are_identical():
    repo_root = Path(__file__).resolve().parents[1]
    root_sty = (repo_root / "mol2lewis.sty").read_text(encoding="utf-8")
    package_sty = (repo_root / "mol2lewis" / "mol2lewis.sty").read_text(encoding="utf-8")

    assert root_sty == package_sty
