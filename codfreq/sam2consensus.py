import json
import cython  # type: ignore
from collections import defaultdict, Counter
from .codfreq_types import (
    NAPos,
    NAChar,
    Profile,
    FragmentConfig,
    SequenceAssemblyConfig,
    NARegionConfig,
    RegionalConsensus
)
from .posnas import get_posnas_in_genome_region

from .filename_helper import name_bamfile


GAP = ord(b'N')
ENCODING = 'UTF-8'


@cython.cfunc
@cython.inline
@cython.returns(dict)
def make_consensus(
    nacons_lookup: dict[tuple[NAPos, int], NAChar],
    region: NARegionConfig
) -> RegionalConsensus:
    """Build consensus sequence for a region from nucleotide counts.

    :param nacons_lookup: Mapping of ``(refpos, insertion_index)`` to the
        nucleotide's ordinal value.
    :type nacons_lookup: dict[tuple[NAPos, int], NAChar]
    :param region: Region definition including name and coordinate range.
    :type region: NARegionConfig
    :returns: Consensus record for the region.
    :rtype: RegionalConsensus
    """

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
    """Generate a consensus sequence for a region from a SAM file.

    :param sampath: Path to the SAM/BAM file.
    :type sampath: str
    :param region: Region configuration describing fragment and coordinates.
    :type region: NARegionConfig
    :returns: Consensus nucleotides covering the region.
    :rtype: RegionalConsensus
    """

    nafreqs: defaultdict[
        tuple[NAPos, int],
        Counter[NAChar]
    ] = defaultdict(Counter)

    for _, posnas in get_posnas_in_genome_region(
        sampath,
        ref_name=region['fromFragment'],
        ref_start=region['refStart'],
        ref_end=region['refEnd']
    ):
        for refpos, idx, na, _ in posnas:
            nafreqs[(refpos, idx)][na] += 1

    nacons_with_count_lookup: dict[tuple[NAPos, int], tuple[NAChar, int]] = {
        pos: nas.most_common(1)[0]
        for pos, nas in nafreqs.items()
    }
    nacons_lookup: dict[tuple[NAPos, int], NAChar] = {
        (pos, idx): na
        for (pos, idx), (na, count) in nacons_with_count_lookup.items()
        if idx == 0 or
        # insertion should only be kept when it's at least as 50% common as the
        # prior nucleotide
        count * 2 > nacons_with_count_lookup.get((pos, 0), ('.', 0))[1]
    }

    r: RegionalConsensus = make_consensus(nacons_lookup, region)
    return r


@cython.ccall
@cython.returns(cython.void)
def create_untrans_region_consensus(
    seqname: str,
    profile: Profile
) -> None:
    """Write consensus sequences for untranslated regions.

    Fragments lacking a ``fromFragment`` source are scanned for regions in the
    profile's ``sequenceAssemblyConfig``. Consensus strings are written to a
    ``<seqname>.untrans.json`` file.

    :param seqname: Base name used to resolve SAM files and the output path.
    :type seqname: str
    :param profile: Profile describing fragments and assembly regions.
    :type profile: Profile
    :rtype: None
    """
    refname: str
    samfile: str
    fragment: FragmentConfig
    region: SequenceAssemblyConfig

    results: list[RegionalConsensus] = []
    for fragment in profile['fragmentConfig']:
        if 'fromFragment' in fragment:
            continue
        refname = fragment['fragmentName']
        samfile = name_bamfile(seqname, refname, is_trimmed=True)
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
    with open(f'{seqname}.untrans.json', 'w') as fp:
        json.dump(results, fp)
