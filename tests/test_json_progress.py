import json
from typing import Any

from codfreq.json_progress import JsonProgress


def test_json_progress_update_and_close(capsys: Any) -> None:
    prog = JsonProgress("task", total=1, ts_interval=0)
    prog.set_description("updated")
    prog.update(1)
    working = json.loads(capsys.readouterr().out.strip())
    assert working["description"] == "updated"
    prog.close()
    done = json.loads(capsys.readouterr().out.strip())
    assert done["status"] == "done"


def test_json_progress_no_output_before_interval(capsys: Any) -> None:
    prog = JsonProgress("task", total=1, ts_interval=10**9)
    prog.update(1)
    assert capsys.readouterr().out == ""
