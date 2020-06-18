import sys
import time
import json
import click


class JsonProgress:

    def __init__(self, description, total, ts_interval=1000, **extras):
        self.description = description
        self.total = total
        self.count = 0
        self.ts_interval = ts_interval
        self.prev_ts = int(time.time() * 1000)
        self.extras = extras

    def update(self, count):
        self.count += count
        now = int(time.time() * 1000)
        if now - self.prev_ts >= self.ts_interval:
            self.prev_ts = now
            click.echo(json.dumps({
                'op': 'progress',
                'status': 'working',
                'description': self.description,
                'count': self.count,
                'total': self.total,
                'ts': now,
                **self.extras
            }))
            sys.stdout.flush()

    def close(self):
        now = int(time.time() * 1000)
        click.echo(json.dumps({
            'op': 'progress',
            'status': 'done',
            'description': self.description,
            'count': self.count,
            'total': self.total,
            'ts': now,
            **self.extras
        }))
