import sys
import time
import json
import click  # type: ignore

from typing import Dict, Any


class JsonProgress:
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
        self.count += count
        now: int = int(time.time() * 1000)
        if now - self.prev_ts >= self.ts_interval:
            self.prev_ts = now
            click.echo(json.dumps({
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
        now: int = int(time.time() * 1000)
        click.echo(json.dumps({
            'op': self.op,
            'status': 'done',
            'description': self.description,
            'count': self.count,
            'total': self.total,
            'ts': now,
            **self.extras
        }))
