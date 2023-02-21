import os
import json
from tqdm import tqdm  # type: ignore
from typing import List, Tuple, Optional, Dict, Iterable, Union, TextIO

from .codfreq_types import (
    DerivedFragmentConfig, MainFragmentConfig,
    Profile,
    FragmentConfig,
    SequenceAssemblyConfig,
    RegionalConsensus
)
from .posnas import PosNA, join_posnas
from .json_progress import JsonProgress
from .segfreq import SegFreq, DEFAULT_CONSENSUS_LEVEL
from .filename_helper import name_segfreq, name_consensus
from .profile import get_ref_fragments
from .fastareader import write_multi_alignment


PATTERNS_HEADER: List[str] = [
    'gene', 'pos', 'nuc',
    'count', 'percent'
]
ENCODING: str = 'UTF-8-sig'


def iter_reference(
    ref: MainFragmentConfig,
    fragments: List[DerivedFragmentConfig]
) -> Iterable[Tuple[str, float, Tuple[Optional[PosNA], ...]]]:
    for fragment in fragments:
        if 'consensus' not in fragment['outputs']:
            continue
        if fragment['geneName'] is None:
            continue
        for level in fragment['outputOptions'].get(
                'consensusLevels', [DEFAULT_CONSENSUS_LEVEL]):
            seq: List[Optional[PosNA]] = []
            for pos_start, pos_end in fragment['refRanges']:
                seq.extend([
                    PosNA(pos_start + idx, 0, ord(na))
                    for idx, na in enumerate(
                        ref['refSequence'][pos_start - 1:pos_end]
                    )
                ])
            yield fragment['geneName'], level, tuple(seq)


def iter_consensus(
    name: str,
    refname: str,
    fragments: List[DerivedFragmentConfig]
) -> Iterable[Tuple[str, float, Tuple[Optional[PosNA], ...]]]:
    segfreq_file: str = name_segfreq(name, refname)
    with open(segfreq_file, encoding=ENCODING) as fh:
        segfreq = SegFreq.load(fh)
        for fragment in fragments:
            if 'consensus' not in fragment['outputs']:
                continue
            if fragment['geneName'] is None:
                continue
            for level in fragment['outputOptions'].get(
                    'consensusLevels', [DEFAULT_CONSENSUS_LEVEL]):
                seq: List[Optional[PosNA]] = []
                for pos_start, pos_end in fragment['refRanges']:
                    seq.extend(
                        segfreq.get_consensus(pos_start, pos_end, level))
                yield fragment['geneName'], level, tuple(seq)


def save_consensus(
    names: List[str],
    profile: Profile,
    log_format: str = 'text'
) -> None:
    ref_fragments = get_ref_fragments(profile)
    if not any(
        'consensus' in frag['outputs']
        for _, _, fragments in ref_fragments
        for frag in fragments
    ):
        return
    groups: Dict[
        Tuple[str, float],
        List[Tuple[str, Tuple[Optional[PosNA], ...]]]
    ] = {}

    total = len(names)
    pbar: Optional[Union[JsonProgress, tqdm]]
    if log_format == 'json':
        pbar = JsonProgress(
            total=total, description='consensus', op='save-consensus')
    elif log_format == 'text':
        pbar = tqdm(total=total)
        pbar.set_description('Generating consensus')
    for refname, ref, fragments in ref_fragments:
        for gene, level, seq in iter_reference(ref, fragments):
            groups[(gene, level)] = [(refname, seq)]

    for name in names:
        basename = os.path.basename(name)
        for refname, _, fragments in ref_fragments:
            for gene, level, seq in iter_consensus(name, refname, fragments):
                groups[(gene, level)].append((
                    '{}|{}|{:g}%'.format(basename, gene, level * 100),
                    seq
                ))
        if pbar:
            pbar.update(1)
    if pbar:
        pbar.close()

    dirname = os.path.dirname(names[0])
    for (gene, level), seqs in groups.items():
        consensus_file: str = name_consensus(dirname, gene, level)
        with open(consensus_file, 'w', encoding='UTF-8') as fh:
            write_multi_alignment(fh, seqs)


def save_untrans_region_consensus(
    name: str,
    profile: Profile,
    log_format: str = 'text'
) -> None:
    refname: str
    samfile: str
    segfreq_file: str
    fragment: FragmentConfig
    region: SequenceAssemblyConfig

    fh: Optional[TextIO] = None
    results: List[RegionalConsensus] = []

    for fragment in profile['fragmentConfig']:
        if 'fromFragment' in fragment:
            continue
        refname = fragment['fragmentName']

        segfreq_file = name_segfreq(name, refname)
        try:
            fh = open(segfreq_file, encoding=ENCODING)
            segfreq = SegFreq.load(fh)

            for region in profile['sequenceAssemblyConfig']:
                if region.get('fromFragment') != refname:
                    continue
                if 'name' not in region or region['name'] is None:
                    continue
                if 'fromFragment' not in region or \
                        region['fromFragment'] is None:
                    continue
                if 'refStart' not in region or region['refStart'] is None:
                    continue
                if 'refEnd' not in region or region['refEnd'] is None:
                    continue
                cons: Tuple[Optional[PosNA]] = segfreq.get_consensus(
                    region['refStart'], region['refEnd'], level=1.)

                results.append(join_posnas(cons))
        finally:
            if fh is not None:
                fh.close()
        with open('{}.untrans.json'.format(name), 'w') as fp:
            json.dump(results, fp)
