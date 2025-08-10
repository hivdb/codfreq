from codfreq.posnas_types import PosNA


def test_posna_alias_structure() -> None:
    posna: PosNA = (1, 0, ord("A"), 30)
    assert posna == (1, 0, ord("A"), 30)
