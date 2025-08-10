from collections import Counter
from types import SimpleNamespace
from typing import Any, Iterator, Literal, cast

import codfreq.sam2codfreq as s2c
from codfreq.codfreq_types import (
    DerivedFragmentConfig,
    MainFragmentConfig,
    Profile,
)


def test_build_fragment_intervals_and_get_ref_ranges() -> None:
    fragments = [
        {"refRanges": [(1, 3), (5, 7)], "fragmentName": "fragA"},
        {"refRanges": [(8, 9)], "fragmentName": "fragB"},
    ]
    assert s2c.build_fragment_intervals(fragments) == [
        ([(1, 3), (5, 7)], "fragA"),
        ([(8, 9)], "fragB"),
    ]

    assert s2c.get_ref_ranges({"refRanges": [(1, 3)]}) == [(1, 3)]
    assert s2c.get_ref_ranges({"refStart": 4, "refEnd": 6}) == [(4, 6)]


def test_get_ref_fragments() -> None:
    profile = {
        "fragmentConfig": [
            {"fragmentName": "refA", "refSequence": "AAA"},
            {
                "fragmentName": "fragA",
                "fromFragment": "refA",
                "geneName": "geneX",
                "refRanges": [(1, 3)],
                "codonAlignment": [{"relRefStart": 1, "windowSize": 3}],
            },
        ]
    }
    refs, lookup = s2c.get_ref_fragments(profile)
    assert refs == [
        (
            "refA",
            {"fragmentName": "refA", "refSequence": "AAA"},
            [
                {
                    "fragmentName": "fragA",
                    "fromFragment": "refA",
                    "geneName": "geneX",
                    "refRanges": [(1, 3)],
                    "codonAlignment": [{"relRefStart": 1, "windowSize": 3}],
                }
            ],
        )
    ]
    assert lookup == {"fragA": [("geneX", 0)]}


def test_to_codon_counter_by_fragpos_and_get_codonfreq() -> None:
    codon_counter = Counter({("fragA", 1, "AAA"): 2, ("fragA", 1, "CCC"): 1})
    fragpos = s2c.to_codon_counter_by_fragpos(codon_counter)
    assert fragpos == {("fragA", 1): Counter({"AAA": 2, "CCC": 1})}

    qualities = {("fragA", 1): Counter({"AAA": 50, "CCC": 20})}
    lookup = {"fragA": [("geneX", 0)]}
    rows = s2c.get_codonfreq(fragpos, qualities, lookup)
    assert rows == [
        {
            "gene": "geneX",
            "position": 1,
            "total": 3,
            "codon": "AAA",
            "count": 2,
            "total_quality_score": 50,
        },
        {
            "gene": "geneX",
            "position": 1,
            "total": 3,
            "codon": "CCC",
            "count": 1,
            "total_quality_score": 20,
        },
    ]


def test_sam2codfreq_between(monkeypatch: Any) -> None:
    def stub_iter_poscodons(
        *args: Any, **kwargs: Any
    ) -> Iterator[tuple[None, list[tuple[str, int, str, int]]]]:
        yield None, [("fragA", 1, "AAA", 30)]
        yield None, [("fragA", 1, "AAA", 20), ("fragA", 2, "CCC", 10)]

    monkeypatch.setattr(s2c, "iter_poscodons", stub_iter_poscodons)
    stat, qual, num_row = s2c.sam2codfreq_between(
        "file.bam", 0, 10, [([(1, 3)], "fragA")]
    )
    assert stat == Counter({("fragA", 1, "AAA"): 2, ("fragA", 2, "CCC"): 1})
    assert qual == Counter({("fragA", 1, "AAA"): 50, ("fragA", 2, "CCC"): 10})
    assert num_row == 2


def test_sam2codfreq(monkeypatch: Any) -> None:
    class DummyAlignmentFile:
        mapped: int

        def __init__(self, *args: Any, **kwargs: Any) -> None:
            self.mapped = 3

        def __enter__(self) -> "DummyAlignmentFile":
            return self

        def __exit__(self, *exc: Any) -> Literal[False]:
            return False

    class DummyTQDM:
        total: int

        def __init__(self, total: int) -> None:
            self.total = total

        def set_description(self, desc: str) -> None:
            self.desc = desc

        def update(self, n: int) -> None:
            self.total += n

        def close(self) -> None:
            pass

    monkeypatch.setattr(
        s2c,
        "pysam",
        SimpleNamespace(AlignmentFile=DummyAlignmentFile),
    )
    monkeypatch.setattr(s2c, "tqdm", DummyTQDM)
    monkeypatch.setattr(
        s2c,
        "chunked_samfile",
        lambda *args, **kwargs: [(0, 1)],
    )

    class DummyExecutor:
        def __init__(self, workers: int) -> None:
            self.workers = workers

        def __enter__(self) -> "DummyExecutor":
            return self

        def __exit__(self, *exc: Any) -> Literal[False]:
            return False

        def map(self, fn: Any, *iterables: Any) -> list[Any]:
            return [fn(*args) for args in zip(*iterables)]

    monkeypatch.setattr(s2c, "ProcessPoolExecutor", DummyExecutor)

    def stub_between(*args: Any, **kwargs: Any) -> tuple[Counter, Counter, int]:
        return (
            Counter({("fragA", 1, "AAA"): 2}),
            Counter({("fragA", 1, "AAA"): 40}),
            2,
        )

    monkeypatch.setattr(s2c, "sam2codfreq_between", stub_between)
    monkeypatch.setattr(
        s2c,
        "codonalign_consensus",
        lambda stat, qual, ref, frags: (stat, qual),
    )

    ref = cast(MainFragmentConfig, {"fragmentName": "refA", "refSequence": "AAA"})
    fragments = [
        cast(
            DerivedFragmentConfig,
            {
                "fragmentName": "fragA",
                "fromFragment": "refA",
                "refRanges": [(1, 3)],
            },
        )
    ]
    stat_by_fragpos, qual_by_fragpos = s2c.sam2codfreq(
        "file.bam", ref, fragments, workers=1
    )
    assert stat_by_fragpos == {("fragA", 1): Counter({"AAA": 2})}
    assert qual_by_fragpos == {("fragA", 1): Counter({"AAA": 40})}


def test_sam2codfreq_all(monkeypatch: Any) -> None:
    def stub_sam2codfreq(
        *args: Any, **kwargs: Any
    ) -> tuple[dict[tuple[str, int], Counter[str]], dict[tuple[str, int], Counter[str]]]:
        return (
            {("fragA", 1): Counter({"AAA": 1})},
            {("fragA", 1): Counter({"AAA": 30})},
        )

    monkeypatch.setattr(s2c, "sam2codfreq", stub_sam2codfreq)
    monkeypatch.setattr(
        s2c, "name_bamfile", lambda name, refname, is_trimmed=True: "file.bam"
    )

    profile = cast(
        Profile,
        {
            "fragmentConfig": [
                {"fragmentName": "refA", "refSequence": "AAA"},
                {
                    "fragmentName": "fragA",
                    "fromFragment": "refA",
                    "geneName": "geneX",
                    "refRanges": [(1, 3)],
                },
            ]
        },
    )

    rows = s2c.sam2codfreq_all("sample", (None, None), profile, workers=1)
    assert rows == [
        {
            "gene": "geneX",
            "position": 1,
            "total": 1,
            "codon": "AAA",
            "count": 1,
            "total_quality_score": 30,
        }
    ]
