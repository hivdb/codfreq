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
from .filename_helper import name_segfreq, name_codfreq
from .profile import get_ref_fragments


CODFREQ_HEADER: List[str] = [
    'gene', 'position',
    'total', 'codon', 'count'
]
ENCODING: str = 'UTF-8-sig'


def iter_codfreq(
    name: str,
    refname: str,
    fragments: List[DerivedFragmentConfig]
) -> Iterable[Dict[str, Any]]:
    segfreq_file: str = name_segfreq(name, refname)
    with open(segfreq_file, encoding=ENCODING) as fh:
        segfreq = SegFreq.load(fh)
        for fragment in fragments:
            if 'codfreq' not in fragment['outputs']:
                continue
            offset = 1
            positions = chain(*[
                range(pos_start, pos_end + 1)
                for pos_start, pos_end in fragment['refRanges']
            ])
            for pos in chunked(positions, 3):
                codfreq = segfreq.get_frequency(pos, 3)
                for codon, count in codfreq.items():
                    yield {
                        'gene': fragment['geneName'],
                        'position': offset,
                        'total': sum(codfreq.values()),
                        'codon': codon.decode('ASCII'),
                        'count': count
                    }
                offset += 1


def save_codfreq(
    name: str,
    fnpair: Tuple[Optional[FASTQFileName], ...],
    profile: Profile,
    log_format: str = 'text'
) -> None:
    ref_fragments = get_ref_fragments(profile)
    if not any(
        'codfreq' in frag['outputs']
        for _, _, fragments in ref_fragments
        for frag in fragments
    ):
        return

    codfreq_file: str = name_codfreq(name)
    with open(codfreq_file, 'w', encoding=ENCODING) as fh:
        writer = csv.DictWriter(fh, CODFREQ_HEADER)
        writer.writeheader()
        for refname, _, fragments in ref_fragments:
            for row in iter_codfreq(
                name, refname, fragments
            ):
                writer.writerow(row)

    if log_format == 'text':
        click.echo(
            'Save codfreq to {}'
            .format(codfreq_file)
        )
    else:
        click.echo(json.dumps({
            'op': 'save-codfreq',
            'status': 'done',
            'query': name,
            'target': codfreq_file
        }))
