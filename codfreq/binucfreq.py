import csv
import json
import click
from typing import List, Tuple, Optional, Dict, Any, Iterable
from itertools import chain
from more_itertools import chunked


from .codfreq_types import (
    FASTQFileName, Profile, DerivedFragmentConfig
)
from .segfreq import SegFreq
from .filename_helper import name_segfreq, name_binucfreq
from .profile import get_ref_fragments


BINUCFREQ_HEADER: List[str] = [
    'gene', 'na_position',
    'total', 'binuc', 'count'
]
ENCODING: str = 'UTF-8-sig'


def iter_binucfreq(
    name: str,
    refname: str,
    fragments: List[DerivedFragmentConfig]
) -> Iterable[Dict[str, Any]]:
    segfreq_file: str = name_segfreq(name, refname)
    with open(segfreq_file, encoding=ENCODING) as fh:
        segfreq = SegFreq.load(fh)
        for fragment in fragments:
            if 'binucfreq' not in fragment['outputs']:
                continue
            offset = 1
            positions = chain(*[
                range(pos_start, pos_end + 1)
                for pos_start, pos_end in fragment['refRanges']
            ])
            for pos in chunked(positions, 2):
                binucfreq = segfreq.get_frequency(pos, 2)
                for binuc, count in binucfreq.items():
                    yield {
                        'gene': fragment['geneName'],
                        'na_position': offset,
                        'total': sum(binucfreq.values()),
                        'binuc': binuc.decode('ASCII'),
                        'count': count
                    }
                offset += 2


def save_binucfreq(
    name: str,
    fnpair: Tuple[Optional[FASTQFileName], ...],
    profile: Profile,
    log_format: str = 'text'
) -> None:
    ref_fragments = get_ref_fragments(profile)
    if not any(
        'binucfreq' in frag['outputs']
        for _, _, fragments in ref_fragments
        for frag in fragments
    ):
        return

    binucfreq_file: str = name_binucfreq(name)
    with open(binucfreq_file, 'w', encoding=ENCODING) as fh:
        writer = csv.DictWriter(fh, BINUCFREQ_HEADER)
        writer.writeheader()
        for refname, _, fragments in ref_fragments:
            for row in iter_binucfreq(
                name, refname, fragments
            ):
                writer.writerow(row)

    if log_format == 'text':
        click.echo(
            'Save binucfreq to {}'
            .format(binucfreq_file)
        )
    else:
        click.echo(json.dumps({
            'op': 'save-binucfreq',
            'status': 'done',
            'query': name,
            'target': binucfreq_file
        }))
