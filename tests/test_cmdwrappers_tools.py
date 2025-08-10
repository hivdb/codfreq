from pathlib import Path
from unittest.mock import MagicMock, patch, call, ANY

import pytest

from codfreq.cmdwrappers import cutadapt, fastp, ivar, minimap2, pigz, samtools


def test_cutadapt_load_config(tmp_path: Path) -> None:
    config = tmp_path / "cutadapt.json"
    config.write_text('{"error_rate": 0.1}', encoding="utf-8")
    adapter3 = tmp_path / "a3.fa"
    adapter3.write_text("", encoding="utf-8")
    adapter5 = tmp_path / "a5.fa"
    adapter5.write_text("", encoding="utf-8")
    adapter53 = tmp_path / "a53.fa"
    adapter53.write_text("", encoding="utf-8")

    loaded = cutadapt.load_config(
        str(config), str(adapter3), str(adapter5), str(adapter53)
    )
    assert loaded == {
        "error_rate": 0.1,
        "adapter3": f"file:{adapter3}",
        "adapter5": f"file:{adapter5}",
        "adapter53": f"file:{adapter53}",
    }
    assert (
        cutadapt.load_config(
            "missing.json",
            str(adapter3),
            str(adapter5),
            str(adapter53),
        )
        is None
    )


def test_cutadapt_runs_and_logs(tmp_path: Path) -> None:
    fake_proc = MagicMock(returncode=0)
    fake_proc.communicate.return_value = ("out", "err")
    out_path = tmp_path / "out.fq"
    with patch.object(cutadapt, "Popen", return_value=fake_proc) as popen_mock:
        cutadapt.cutadapt(
            "in.fq",
            str(out_path),
            adapter3="A3",
            adapter53="B",
            error_rate=0.2,
            no_indels=True,
            times=2,
            min_overlap=3,
        )
    expected_cmd = [
        "cutadapt",
        "-j",
        "0",
        "-a",
        "A3",
        "-b",
        "B",
        "-e",
        "0.2",
        "--no-indels",
        "-n",
        "2",
        "-O",
        "3",
        "-o",
        str(out_path),
        "in.fq",
    ]
    popen_mock.assert_called_once_with(
        expected_cmd, stdout=cutadapt.PIPE, stderr=cutadapt.PIPE, encoding="U8"
    )
    log_file = out_path.with_suffix(".fq.cutadapt.log")
    assert log_file.exists()
    assert "cutadapt" in log_file.read_text()


def test_fastp_load_config(tmp_path: Path) -> None:
    config = tmp_path / "fastp.json"
    config.write_text('{"qualified_quality_phred": 15}', encoding="utf-8")
    loaded = fastp.load_config(str(config))
    assert loaded["qualified_quality_phred"] == 15
    assert fastp.load_config(str(tmp_path / "missing.json")) == {}


def test_fastp_paired_parameters(tmp_path: Path) -> None:
    fake_proc = MagicMock(returncode=0)
    fake_proc.communicate.return_value = ("out", "err")
    merged = tmp_path / "merged.fq"
    with patch.object(fastp, "Popen", return_value=fake_proc) as popen_mock:
        fastp.fastp(
            "r1.fq",
            "r2.fq",
            str(merged),
            include_unmerged=True,
            qualified_quality_phred=15,
            unqualified_percent_limit=30,
            n_base_limit=5,
            average_qual=20,
            length_required=50,
            length_limit=100,
            adapter_sequence="AAA",
            adapter_sequence_r2="TTT",
        )
    expected = [
        "fastp",
        "-w",
        fastp.THREADS,
        "-i",
        "r1.fq",
        "-I",
        "r2.fq",
        "-m",
        "--merged_out",
        str(merged),
        "--include_unmerged",
        "-q",
        "15",
        "-u",
        "30",
        "-n",
        "5",
        "-e",
        "20",
        "-l",
        "50",
        "--length_limit",
        "100",
        "-a",
        "AAA",
        "--adapter_sequence_r2",
        "TTT",
    ]
    popen_mock.assert_called_once_with(
        expected, stdout=fastp.PIPE, stderr=fastp.PIPE, encoding="U8"
    )
    assert (tmp_path / "merged.fastp.log").exists()


def test_fastp_single_end_disable(tmp_path: Path) -> None:
    fake_proc = MagicMock(returncode=0)
    fake_proc.communicate.return_value = ("out", "err")
    out = tmp_path / "out.fq"
    with patch.object(fastp, "Popen", return_value=fake_proc) as popen_mock:
        fastp.fastp(
            "r1.fq",
            None,
            str(out),
            disable_quality_filtering=True,
            disable_length_filtering=True,
            disable_adapter_trimming=True,
        )
    expected = [
        "fastp",
        "-w",
        fastp.THREADS,
        "-i",
        "r1.fq",
        "-o",
        str(out),
        "-Q",
        "-L",
        "-A",
    ]
    popen_mock.assert_called_once_with(
        expected, stdout=fastp.PIPE, stderr=fastp.PIPE, encoding="U8"
    )
    assert (tmp_path / "out.fastp.log").exists()


def test_ivar_load_trim_config(tmp_path: Path) -> None:
    config = tmp_path / "trim.json"
    config.write_text('{"min_length": 50}', encoding="utf-8")
    primers = tmp_path / "primers.bed"
    primers.write_text("", encoding="utf-8")
    loaded = ivar.load_trim_config(str(config), str(primers))
    assert loaded == {"min_length": 50, "primers_bed": str(primers)}
    assert ivar.load_trim_config("missing.json", str(primers)) is None


def test_ivar_trim_runs(tmp_path: Path) -> None:
    input_bam = tmp_path / "in.bam"
    input_bam.write_text("", encoding="utf-8")
    output_bam = tmp_path / "out.bam"
    fake_proc = MagicMock(returncode=0)
    fake_proc.communicate.return_value = ("out", "err")
    with (
        patch.object(ivar, "Popen", return_value=fake_proc) as popen_mock,
        patch.object(
            ivar,
            "execute",
            side_effect=[("sortout", "sorterr"), ("idxout", "idxerr")],
        ) as exec_mock,
    ):
        ivar.trim(
            str(input_bam),
            str(output_bam),
            primers_bed="primers.bed",
            min_length=50,
            include_reads_with_no_primers=True,
        )
    command = popen_mock.call_args[0][0]
    assert command[0:4] == ["ivar", "trim", "-i", str(input_bam)]
    assert "-b" in command and "primers.bed" in command
    assert "-m" in command and "50" in command
    assert "-e" in command
    exec_mock.assert_has_calls(
        [
            call(
                [
                    "samtools",
                    "sort",
                    "-@",
                    ivar.THREADS,
                    "-O",
                    "bam",
                    "-o",
                    str(output_bam),
                    ANY,
                ]
            ),
            call(["samtools", "index", "-@", ivar.THREADS, str(output_bam)]),
        ]
    )
    assert (tmp_path / "out.bam.ivar.log").exists()


def test_minimap2_refinit_registration() -> None:
    from codfreq.cmdwrappers.base import REFINIT_FUNCTIONS

    assert "minimap2" in REFINIT_FUNCTIONS
    assert minimap2.minimap2_refinit("ref.fa") is None


def test_minimap2_align_runs(tmp_path: Path) -> None:
    bam = tmp_path / "out.bam"
    proc_minimap2 = MagicMock(returncode=0)
    proc_minimap2.stdout = MagicMock()
    proc_minimap2.stderr = MagicMock()
    proc_minimap2.stderr.read.return_value = "minimap2 err"
    proc_sam2bam = MagicMock(returncode=0)
    proc_sam2bam.communicate.return_value = ("out", "err")
    with (
        patch.object(
            minimap2, "Popen", side_effect=[proc_minimap2, proc_sam2bam]
        ) as popen_mock,
        patch.object(
            minimap2, "execute", return_value=("idxout", "idxerr")
        ) as exec_mock,
    ):
        result = minimap2.minimap2_align("ref.fa", "r1.fq", "r2.fq", str(bam))
    expected_align = [
        "minimap2",
        *minimap2.MINIMAP2_ARGS,
        "-a",
        "ref.fa",
        "r1.fq",
        "r2.fq",
    ]
    expected_sort = [
        "samtools",
        "sort",
        "-@",
        minimap2.THREADS,
        "-O",
        "bam",
        "-o",
        str(bam),
    ]
    popen_mock.assert_has_calls(
        [
            call(
                expected_align,
                stdout=minimap2.PIPE,
                stderr=minimap2.PIPE,
                encoding="U8",
            ),
            call(
                expected_sort,
                stdin=proc_minimap2.stdout,
                stdout=minimap2.PIPE,
                stderr=minimap2.PIPE,
                encoding="U8",
            ),
        ]
    )
    exec_mock.assert_called_once_with(
        ["samtools", "index", "-@", minimap2.THREADS, str(bam)]
    )
    assert (tmp_path / "out.bam.minimap2.log").exists()
    assert result == {"overall_rate": -1.}


def test_pigz_compress_and_error() -> None:
    good_proc = MagicMock(returncode=0)
    good_proc.communicate.return_value = (b"data", b"")
    with patch.object(pigz, "Popen", return_value=good_proc) as popen_mock:
        out = pigz.compress(b"in", compresslevel=5, mtime=1.0)
    popen_mock.assert_called_once_with(
        ["pigz", "-5", "-c", "-M", "1.0"],
        stdin=pigz.PIPE,
        stdout=pigz.PIPE,
        stderr=pigz.PIPE,
    )
    assert out == b"data"
    err_proc = MagicMock(returncode=0)
    err_proc.communicate.return_value = (b"", b"fail")
    with patch.object(pigz, "Popen", return_value=err_proc):
        with pytest.raises(RuntimeError):
            pigz.compress(b"in")


def test_pigz_decompress() -> None:
    fake_proc = MagicMock(returncode=0)
    fake_proc.communicate.return_value = (b"out", b"")
    with patch.object(pigz, "Popen", return_value=fake_proc) as popen_mock:
        out = pigz.decompress(b"in")
    popen_mock.assert_called_once_with(
        ["pigz", "-d", "-c"],
        stdin=pigz.PIPE,
        stdout=pigz.PIPE,
        stderr=pigz.PIPE,
    )
    assert out == b"out"


def test_samtools_stats(tmp_path: Path) -> None:
    sam = tmp_path / "input.sam"
    sam.write_text("", encoding="utf-8")
    with patch.object(
        samtools, "execute", return_value=("logs", "")
    ) as exec_mock:
        samtools.stats(str(sam))
    exec_mock.assert_called_once_with(["samtools", "stats", str(sam)])
    stats_file = tmp_path / "input.stats.txt"
    assert stats_file.read_text() == "logs"
