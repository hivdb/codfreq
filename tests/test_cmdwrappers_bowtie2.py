"""Tests for the Bowtie2 command wrappers."""

from unittest.mock import patch, mock_open
from pathlib import Path

from codfreq.cmdwrappers import bowtie2


def test_bowtie2_refinit_executes_when_missing() -> None:
    """Run ``bowtie2-build`` when index files are missing."""
    with (
        patch(
            "codfreq.cmdwrappers.bowtie2.os.path.isfile",
            return_value=False,
        ),
        patch("codfreq.cmdwrappers.bowtie2.execute") as exec_mock,
    ):
        bowtie2.bowtie2_refinit("ref.fa")
    exec_mock.assert_called_once_with(['bowtie2-build', 'ref.fa', 'ref'])


def test_bowtie2_refinit_skips_when_present() -> None:
    """Existing index files skip build step."""
    with (
        patch(
            "codfreq.cmdwrappers.bowtie2.os.path.isfile",
            return_value=True,
        ),
        patch("codfreq.cmdwrappers.bowtie2.execute") as exec_mock,
    ):
        bowtie2.bowtie2_refinit("ref.fa")
    exec_mock.assert_not_called()


def test_bowtie2_align_parses_rate_and_logs(tmp_path: Path) -> None:
    """Alignment writes a log file and extracts overall rate."""
    logs = "reads\n95.5% overall alignment rate\n"
    with (
        patch(
            "codfreq.cmdwrappers.bowtie2.execute", return_value=(logs, "")
        ) as exec_mock,
        patch(
            "codfreq.cmdwrappers.bowtie2.open", mock_open(), create=True
        ) as m_open,
    ):
        result = bowtie2.bowtie2_align(
            "ref.fa", "r1.fq", "r2.fq", str(tmp_path / "out.sam"))
    exec_mock.assert_called_once()
    m_open.assert_called_once_with(str(tmp_path / "out.log"), 'w')
    m_open().write.assert_called_once_with(logs)
    assert result == {'overall_rate': 95.5}
    cmd = exec_mock.call_args.args[0]
    assert cmd[:3] == ['bowtie2', '--local', '--threads']
    assert cmd[-6:] == [
        '-S',
        str(tmp_path / 'out.sam'),
        '-U',
        'r1.fq',
        '-U',
        'r2.fq',
    ]
