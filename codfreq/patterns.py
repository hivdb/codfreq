import os
import json
import math
import click
from typing import List, Tuple, Optional, Iterable


from .codfreq_types import (
    FASTQFileName, Profile, MainFragmentConfig, DerivedFragmentConfig
)
from .segfreq import SegFreq, DEFAULT_TOP_N_SEEDS
from .filename_helper import name_segfreq, name_patterns
from .profile import get_ref_fragments
from .posnas import PosNA
from .fastareader import write_multi_alignment


PATTERNS_HEADER: List[str] = [
    'gene', 'pos', 'nuc',
    'count', 'percent'
]
ENCODING: str = 'UTF-8-sig'


def iter_patterns(
    name: str,
    refname: str,
    ref: MainFragmentConfig,
    fragment: DerivedFragmentConfig
) -> Iterable[Tuple[str, Tuple[Optional[PosNA], ...]]]:
    segfreq_file: str = name_segfreq(name, refname)
    basename = os.path.basename(name)
    with open(segfreq_file, encoding=ENCODING) as fh:
        segfreq = SegFreq.load(fh)
        for pos_start, pos_end in fragment['refRanges']:
            yield (
                refname,
                tuple([
                    PosNA(pos_start + idx, 0, ord(na))
                    for idx, na in enumerate(
                        ref['refSequence'][pos_start - 1:pos_end]
                    )
                ])
            )
            patmap = segfreq.get_patterns(
                pos_start, pos_end,
                top_n_seeds=(
                    fragment['outputOptions']
                    .get('patternsTopNSeeds', DEFAULT_TOP_N_SEEDS)
                )
            )
            digits = int(math.log10(len(patmap) or 1)) + 1
            for idx, (pattern, (count, pcnt)) in enumerate(patmap.items()):
                yield (
                    '{{}}.{{:0{}d}}|count={{}}|pcnt={{:g}}%'
                    .format(digits)
                    .format(basename, idx + 1, count, pcnt * 100),
                    pattern
                )


def save_patterns(
    name: str,
    fnpair: Tuple[Optional[FASTQFileName], ...],
    profile: Profile,
    log_format: str = 'text'
) -> None:
    ref_fragments = get_ref_fragments(profile)

    for refname, ref, fragments in ref_fragments:
        for fragment in fragments:
            if 'patterns' not in fragment['outputs']:
                continue
            patterns_file: str = name_patterns(name, fragment['fragmentName'])
            patterns = list(iter_patterns(name, refname, ref, fragment))
            with open(patterns_file, 'w', encoding='UTF-8') as fh:
                write_multi_alignment(fh, patterns)

            if log_format == 'text':
                click.echo(
                    'Save {} patterns to {}'
                    .format(len(patterns) - 1, patterns_file)
                )
            else:
                click.echo(json.dumps({
                    'op': 'save-patterns',
                    'total': len(patterns) - 1,
                    'status': 'done',
                    'query': name,
                    'target': patterns_file
                }))
