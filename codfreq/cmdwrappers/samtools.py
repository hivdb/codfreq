import os
from .base import execute


def stats(sam: str) -> None:
    logs: str
    logs, _ = execute(['samtools', 'stats', sam])
    with open(os.path.splitext(sam)[0] + '.stats.txt', 'w') as fp:
        fp.write(logs)
