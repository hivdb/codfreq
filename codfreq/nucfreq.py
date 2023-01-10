import csv
import json
import click
from typing import List, Tuple, Optional, Dict, Any, Iterable


from .codfreq_types import (
    FASTQFileName, Profile, DerivedFragmentConfig
)
from .segfreq import SegFreq
from .filename_helper import name_segfreq, name_nucfreq
from .profile import get_ref_fragments


NUCFREQ_HEADER: List[str] = [
    'gene', 'position',
    'total', 'nuc', 'count'
]
ENCODING: str = 'UTF-8-sig'


def iter_nucfreq(
    name: str,
    refname: str,
    fragments: List[DerivedFragmentConfig]
) -> Iterable[Dict[str, Any]]:
    segfreq_file: str = name_segfreq(name, refname)
    with open(segfreq_file, encoding=ENCODING) as fh:
        segfreq = SegFreq.load(fh)
        for fragment in fragments:
            offset = 1
            for pos_start, pos_end in fragment['refRanges']:
                for pos in range(pos_start, pos_end + 1):
                    nucfreq = segfreq.get_pos_nas(pos)
                    for nuc, count in nucfreq.items():
                        yield {
                            'gene': fragment['geneName'],
                            'position': pos - pos_start + offset,
                            'total': sum(nucfreq.values()),
                            'nuc': nuc.decode('ASCII'),
                            'count': count
                        }
                offset += pos_end - pos_start + 1


def save_nucfreq(
    name: str,
    fnpair: Tuple[Optional[FASTQFileName], ...],
    profile: Profile,
    log_format: str = 'text'
) -> None:
    ref_fragments = get_ref_fragments(profile)
    if not any(
        frag['frequencyType'] == 'nucleotide'
        for _, _, fragments in ref_fragments
        for frag in fragments
    ):
        return

    nucfreq_file: str = name_nucfreq(name)
    with open(nucfreq_file, 'w', encoding=ENCODING) as fh:
        writer = csv.DictWriter(fh, NUCFREQ_HEADER)
        writer.writeheader()
        for refname, _, fragments in ref_fragments:
            for row in iter_nucfreq(
                name, refname, fragments
            ):
                writer.writerow(row)

    if log_format == 'text':
        click.echo(
            'Save nucfreq to {}'
            .format(nucfreq_file)
        )
    else:
        click.echo(json.dumps({
            'op': 'save-nucleotide',
            'status': 'done',
            'query': name,
            'target': nucfreq_file
        }))
