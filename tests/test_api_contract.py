from mol2lewis import core


def test_lewis_returns_empty_list_when_generation_fails(monkeypatch):
    monkeypatch.setattr(core, "_generate_chemfig", lambda *args, **kwargs: None)
    monkeypatch.setattr(core.pcp, "get_compounds", lambda *args, **kwargs: [])

    result = core.lewis("CCO")

    assert isinstance(result, list)
    assert result == []


def test_lewis_invalid_selection_returns_empty_list(monkeypatch):
    class DummyCompound:
        smiles = "CCO"

    monkeypatch.setattr(core.pcp, "get_compounds", lambda *args, **kwargs: [DummyCompound()])

    result = core.lewis("C2H6O", selection="unknown")

    assert isinstance(result, list)
    assert result == []
