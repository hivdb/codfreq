import json
import cython  # type: ignore
from collections import defaultdict, Counter

from typing import (
    DefaultDict,
    Counter as tCounter,
    Tuple,
    List,
    Dict
)
from .codfreq_types import (
    NAPos,
    NAChar,
    Profile,
    FragmentConfig,
    SequenceAssemblyConfig,
    NARegionConfig,
    RegionalConsensus
)
from .posnas import get_posnas_in_genome_region, PosNA

from .filename_helper import name_bamfile


GAP = ord(b'N')
ENCODING = 'UTF-8'


@cython.cfunc
@cython.inline
@cython.returns(dict)
def make_consensus(
    nacons_lookup: Dict[Tuple[NAPos, int], NAChar],
    region: NARegionConfig
) -> RegionalConsensus:
    refpos: NAPos
    idx: int

    name: str = region['name']
    refpos_start: NAPos = region['refStart']
    refpos_end: NAPos = region['refEnd']
    consarr: bytearray = bytearray()

    for refpos in range(refpos_start, refpos_end + 1):
        idx = 0
        while True:
            try:
                na = nacons_lookup[(refpos, idx)]
            except KeyError:
                if idx == 0:
                    consarr.append(GAP)
                break
            consarr.append(na)
            idx += 1
    return {
        'name': name,
        'refStart': refpos_start,
        'refEnd': refpos_end,
        'consensus': consarr.decode(ENCODING)
    }


def sam2consensus(
    sampath: str,
    region: NARegionConfig,
) -> RegionalConsensus:

    name: str
    refpos_start: NAPos
    refpos_end: NAPos
    consarr: bytearray
    posnas: List[PosNA]
    refpos: NAPos
    idx: int
    na: NAChar

    nafreqs: DefaultDict[
        Tuple[NAPos, int],
        tCounter[NAChar]
    ] = defaultdict(Counter)

    for _, posnas in get_posnas_in_genome_region(
        sampath,
        ref_name=region['fromFragment'],
        ref_start=region['refStart'],
        ref_end=region['refEnd']
    ):
        for refpos, idx, na, _ in posnas:
            nafreqs[(refpos, idx)][na] += 1

    nacons_lookup: Dict[Tuple[NAPos, int], NAChar] = {
        pos: nas.most_common(1)[0][0]
        for pos, nas in nafreqs.items()
    }

    r: RegionalConsensus = make_consensus(nacons_lookup, region)
    return r


@cython.ccall
@cython.returns(cython.void)
def create_untrans_region_consensus(
    seqname: str,
    profile: Profile
) -> None:

    refname: str
    samfile: str
    fragment: FragmentConfig
    region: SequenceAssemblyConfig

    results: List[RegionalConsensus] = []
    for fragment in profile['fragmentConfig']:
        if 'fromFragment' in fragment:
            continue
        refname = fragment['fragmentName']
        samfile = name_bamfile(seqname, refname)
        for region in profile['sequenceAssemblyConfig']:
            if region.get('fromFragment') != refname:
                continue
            if 'name' not in region or region['name'] is None:
                continue
            if 'fromFragment' not in region or region['fromFragment'] is None:
                continue
            if 'refStart' not in region or region['refStart'] is None:
                continue
            if 'refEnd' not in region or region['refEnd'] is None:
                continue

            results.append(sam2consensus(
                samfile,
                {
                    'name': region['name'],
                    'fromFragment': region['fromFragment'],
                    'refStart': region['refStart'],
                    'refEnd': region['refEnd']
                }
            ))
    with open('{}.untrans.json'.format(seqname), 'w') as fp:
        json.dump(results, fp)
