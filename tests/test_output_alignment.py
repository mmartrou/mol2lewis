import shutil

import pytest

from mol2lewis import lewis


pytestmark = pytest.mark.skipif(
    shutil.which("mol2chemfig") is None,
    reason="mol2chemfig executable is required for topology integration tests",
)


def test_no_atom_numbers_mode_keeps_cycle_closures_without_phantom_for_4932998():
    result = lewis("c1ccccc1", show_carbons=False, aromatic_circles=True)

    assert result
    normal = result[0]["normal"]
    draft = result[0]["draft"]

    assert "\\phantom{" not in normal
    assert "\\phantom{" not in draft
