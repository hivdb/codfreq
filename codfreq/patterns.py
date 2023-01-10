import csv
import json
import click
from typing import List, Tuple, Optional, Dict, Any, Iterable


from .codfreq_types import (
    FASTQFileName, Profile, DerivedFragmentConfig
)
from .segfreq import SegFreq
from .filename_helper import name_segfreq, name_patterns
from .profile import get_ref_fragments


PATTERNS_HEADER: List[str] = [
    'gene', 'pos', 'nuc',
    'count', 'percent'
]
ENCODING: str = 'UTF-8-sig'


def iter_patterns(
    name: str,
    refname: str,
    fragments: List[DerivedFragmentConfig]
) -> Iterable[Dict[str, Any]]:
    segfreq_file: str = name_segfreq(name, refname)
    with open(segfreq_file, encoding=ENCODING) as fh:
        segfreq = SegFreq.load(fh)
        for fragment in fragments:
            if 'patterns' not in fragment['outputs']:
                continue
            offset = 1
            for pos_start, pos_end in fragment['refRanges']:
                # TODO: allow to config top_n_seeds
                patterns = segfreq.get_patterns(
                    pos_start, pos_end, top_n_seeds=100)
                for pattern, (count, pcnt) in patterns.items():
                    pos_text = ':'.join([
                        '{}'.format(node.pos - pos_start + offset)
                        for node in pattern
                    ])
                    nuc_text = bytes([
                        node.na for node in pattern
                    ]).decode('ASCII')

                    yield {
                        'gene': fragment['geneName'],
                        'pos': pos_text,
                        'nuc': nuc_text,
                        'count': count,
                        'percent': pcnt
                    }
                offset += pos_end - pos_start + 1


def save_patterns(
    name: str,
    fnpair: Tuple[Optional[FASTQFileName], ...],
    profile: Profile,
    log_format: str = 'text'
) -> None:
    ref_fragments = get_ref_fragments(profile)
    if not any(
        'patterns' in frag['outputs']
        for _, _, fragments in ref_fragments
        for frag in fragments
    ):
        return

    patterns_file: str = name_patterns(name)
    with open(patterns_file, 'w', encoding=ENCODING) as fh:
        writer = csv.DictWriter(fh, PATTERNS_HEADER)
        writer.writeheader()
        for refname, _, fragments in ref_fragments:
            for row in iter_patterns(
                name, refname, fragments
            ):
                writer.writerow(row)

    if log_format == 'text':
        click.echo(
            'Save patterns to {}'
            .format(patterns_file)
        )
    else:
        click.echo(json.dumps({
            'op': 'save-patterns',
            'status': 'done',
            'query': name,
            'target': patterns_file
        }))
