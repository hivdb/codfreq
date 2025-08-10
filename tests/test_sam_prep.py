import sys
import types
from collections import Counter
from typing import Any, List

pysam_stub = types.ModuleType("pysam")


class AlignmentFile:  # pragma: no cover - minimal stub
    def __init__(self, *args: Any, **kwargs: Any) -> None:  # noqa: D401
        pass

    def __enter__(self) -> "AlignmentFile":  # noqa: D401
        return self

    def __exit__(self, *exc_info: Any) -> None:  # noqa: D401
        return None

    def fetch(self) -> List[Any]:  # noqa: D401
        return []

    def tell(self) -> int:  # noqa: D401
        return 0

    def seek(self, pos: int) -> None:  # noqa: D401
        pass

    def write(self, read: object) -> None:  # noqa: D401
        pass


def _install_stub() -> None:
    pysam_stub.AlignmentFile = AlignmentFile  # type: ignore[attr-defined]
    sys.modules.setdefault("pysam", pysam_stub)


_install_stub()

from codfreq.sam_prep import count_indel_positions, squash_gaps  # noqa: E402


def test_squash_gaps_merges_indels() -> None:
    cig = ((0, 5), (1, 2), (0, 3), (2, 1), (0, 4))
    assert squash_gaps(cig) == ((0, 5), (1, 1), (0, 7))


def test_count_indel_positions() -> None:
    counter: Counter[int] = Counter()
    cig = ((0, 5), (1, 2), (0, 3), (2, 1), (0, 4))
    count_indel_positions(cig, 100, counter)
    assert counter == Counter({105: 1, 108: 1})
