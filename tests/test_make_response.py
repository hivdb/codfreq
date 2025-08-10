import json
from pathlib import Path
from datetime import datetime

from typer.testing import CliRunner

from codfreq.make_response import yield_codfreqs, yield_untrans, app


def test_make_response_helpers_and_cli(tmp_path: Path) -> None:
    cf = tmp_path / "sample.codfreq"
    cf.write_text(
        "gene,position,total,codon,count,total_quality_score\n"
        "S,1,100,AAA,10,1000\n"
        "S,2,100,AA,5,500\n",
        encoding="utf-8",
    )
    ut = tmp_path / "sample.untrans.json"
    ut_data = [{"name": "UTR1", "refStart": 1, "refEnd": 2, "consensus": "AA"}]
    ut.write_text(json.dumps(ut_data), encoding="utf-8")

    gen = yield_codfreqs(tmp_path)
    name, rows = next(gen)
    row = next(rows)
    assert name == "sample"
    assert row["codon"] == "AAA"
    next(gen, None)

    names_untrans = dict(yield_untrans(tmp_path))
    assert names_untrans["sample"] == ut_data

    runner = CliRunner()
    result = runner.invoke(app, [str(tmp_path), "prefix"])
    assert result.exit_code == 0
    response = json.loads((tmp_path / "response.json").read_text())
    datetime.fromisoformat(response["lastUpdatedAt"])
    assert response["codfreqs"][0]["name"] == "sample.codfreq"
    assert response["codfreqs"][0]["untranslatedRegions"] == ut_data
    assert response["codfreqs"][0]["allReads"][0]["position"] == 1
    assert len(response["codfreqs"][0]["allReads"]) == 1
