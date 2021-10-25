import os
from .base import execute


def stats(sam):
    logs, _ = execute(['samtools', 'stats', sam])
    with open(os.path.splitext(sam)[0] + '.stats.txt', 'wb') as fp:
        fp.write(logs)
