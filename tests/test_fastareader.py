"""Tests for :mod:`codfreq.fastareader`."""

from io import StringIO

from codfreq import fastareader


def test_load_reads_sequences() -> None:
    """``fastareader.load`` parses headers and sequences."""

    data = ">seq1\nATcg\n#comment\n>seq2\nGG\n"
    sequences = fastareader.load(StringIO(data))
    assert [s.model_dump() for s in sequences] == [
        {"header": "seq1", "sequence": "ATCG"},
        {"header": "seq2", "sequence": "GG"},
    ]
