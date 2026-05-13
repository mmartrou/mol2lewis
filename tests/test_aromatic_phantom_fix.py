from mol2lewis import lewis


def test_raw_mol2chemfig_topology_is_preserved_for_aromatic_cases():
    result = lewis(4932998, show_carbons=False, aromatic_circles=True)

    assert result
    normal = result[0]["normal"]
    draft = result[0]["draft"]

    assert "\\phantom{" not in normal
    assert "\\phantom{" not in draft
    assert "\\mcfcringle" in normal
    assert "\\mcfcringle" in draft


def test_draft_only_changes_circle_color_not_structure():
    result = lewis(5442468, show_carbons=False, aromatic_circles=True)

    assert result
    normal = result[0]["normal"]
    draft = result[0]["draft"]

    assert normal.count("\\mcfcringle") == draft.count("\\mcfcringle")
    assert "\\color{red}{\\mcfcringle" in draft
