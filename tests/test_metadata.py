from pathlib import Path
import re

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    tomllib = None

from mol2lewis import __version__


def test_package_version_matches_pyproject():
    repo_root = Path(__file__).resolve().parents[1]
    pyproject_text = (repo_root / "pyproject.toml").read_text(encoding="utf-8")

    if tomllib is not None:
        pyproject = tomllib.loads(pyproject_text)
        version = pyproject["project"]["version"]
    else:
        project_section = re.search(r"\[project\](.*?)(?:\n\[|\Z)", pyproject_text, re.S)
        assert project_section is not None
        match = re.search(r'^version\s*=\s*"([^"]+)"', project_section.group(1), re.M)
        assert match is not None
        version = match.group(1)

    assert version == __version__


def test_sty_files_are_identical():
    repo_root = Path(__file__).resolve().parents[1]
    root_sty = (repo_root / "mol2lewis.sty").read_text(encoding="utf-8")
    package_sty = (repo_root / "mol2lewis" / "mol2lewis.sty").read_text(encoding="utf-8")

    assert root_sty == package_sty
