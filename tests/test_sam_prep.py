import sys
import types
from collections import Counter
from typing import Any, Iterator, List

pysam_stub = types.ModuleType("pysam")


class AlignmentFile:  # pragma: no cover - minimal stub
    """Simplified :class:`pysam.AlignmentFile` for testing."""

    reads_in: List[Any] = []
    written: List[Any] = []

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        self.mode = args[1] if len(args) > 1 else kwargs.get("mode", "r")
        self._index = 0

    def __enter__(self) -> "AlignmentFile":
        return self

    def __exit__(self, *exc_info: Any) -> None:
        return None

    def __iter__(self) -> Iterator[Any]:
        if self.mode.startswith("r"):
            for read in AlignmentFile.reads_in[self._index:]:
                yield read
                self._index += 1

    def fetch(self) -> List[Any]:  # noqa: D401
        return []

    def tell(self) -> int:
        return self._index

    def seek(self, pos: int) -> None:
        self._index = pos

    def write(self, read: object) -> None:
        AlignmentFile.written.append(read)


def _install_stub() -> None:
    pysam_stub.AlignmentFile = AlignmentFile  # type: ignore[attr-defined]
    sys.modules["pysam"] = pysam_stub


_install_stub()

from codfreq.sam_prep import (  # noqa: E402
    count_indel_positions,
    squash_gaps,
)


def test_squash_gaps_merges_indels() -> None:
    cig = ((0, 5), (1, 2), (0, 3), (2, 1), (0, 4))
    assert squash_gaps(cig) == ((0, 5), (1, 1), (0, 7))


def test_count_indel_positions() -> None:
    counter: Counter[int] = Counter()
    cig = ((0, 5), (1, 2), (0, 3), (2, 1), (0, 4))
    count_indel_positions(cig, 100, counter)
    assert counter == Counter({105: 1, 108: 1})


def test_squash_gaps_switches_operation() -> None:
    """Insertions longer than nearby deletions flip the operation."""

    cig = ((0, 5), (1, 2), (0, 3), (2, 4), (0, 1))
    assert squash_gaps(cig) == ((0, 5), (2, 2), (0, 4))


def test_prepare_sam_processes_reads() -> None:
    """Mapped reads receive squashed cigars and are written out."""

    class Read:
        def __init__(
            self,
            unmapped: bool,
            cigar: tuple[tuple[int, int], ...],
            ref_start: int,
        ) -> None:
            self.is_unmapped = unmapped
            self.cigartuples = list(cigar)
            self.reference_start = ref_start

    unmapped = Read(True, (), 0)
    mapped = Read(False, ((0, 5), (1, 2), (0, 3), (2, 1), (0, 4)), 100)
    _install_stub()
    import importlib
    import codfreq.sam_prep as sam_prep_module
    spm = importlib.reload(sam_prep_module)

    alignment_file = spm.AlignmentFile  # type: ignore[attr-defined]
    alignment_file.reads_in = [unmapped, mapped]  # type: ignore[attr-defined]
    alignment_file.written = []  # type: ignore[attr-defined]
    spm.prepare_sam("in.sam", "out.sam")
    written = alignment_file.written  # type: ignore[attr-defined]
    assert len(written) == 2
    assert written[1].cigartuples == [  # type: ignore[attr-defined]
        (0, 5),
        (1, 1),
        (0, 7),
    ]
