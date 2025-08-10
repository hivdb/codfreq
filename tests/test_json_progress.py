import json

from codfreq.json_progress import JsonProgress


def test_json_progress_update_and_close(capsys) -> None:
    prog = JsonProgress("task", total=1, ts_interval=0)
    prog.set_description("updated")
    prog.update(1)
    working = json.loads(capsys.readouterr().out.strip())
    assert working["description"] == "updated"
    prog.close()
    done = json.loads(capsys.readouterr().out.strip())
    assert done["status"] == "done"
