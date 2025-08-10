import sys
import time
import json
import rich  # type: ignore[import-not-found]

from typing import Dict, Any


class JsonProgress:
    """Report progress to STDOUT in JSON format."""

    description: str
    total: int
    count: int
    prev_ts: int
    ts_interval: int
    op: str
    extras: Dict[str, Any]

    def __init__(
        self,
        description: str,
        total: int,
        ts_interval: int = 1000,
        op: str = 'progress',
        **extras: Any
    ):
        self.description = description
        self.total = total
        self.count = 0
        self.ts_interval = ts_interval
        self.prev_ts = int(time.time() * 1000)
        self.op = op
        self.extras = extras

    def update(self, count: int) -> None:
        """
        Update progress counters and emit a working status when the
        timestamp interval has elapsed.

        :param count: Number of completed units since the last update.
        :returns: None
        """
        self.count += count
        now: int = int(time.time() * 1000)
        if now - self.prev_ts >= self.ts_interval:
            self.prev_ts = now
            rich.print(json.dumps({
                'op': self.op,
                'status': 'working',
                'description': self.description,
                'count': self.count,
                'total': self.total,
                'ts': now,
                **self.extras
            }))
            sys.stdout.flush()

    def close(self) -> None:
        """
        Emit a final message indicating that the task is complete.

        :returns: None
        """
        now: int = int(time.time() * 1000)
        rich.print(json.dumps({
            'op': self.op,
            'status': 'done',
            'description': self.description,
            'count': self.count,
            'total': self.total,
            'ts': now,
            **self.extras
        }))

    def set_description(self, description: str) -> None:
        """Set a new description for subsequent progress updates.

        :param description: Text describing the current task.
        :returns: None
        """
        self.description = description
