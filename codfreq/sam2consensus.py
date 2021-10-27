import json
from collections import defaultdict, Counter

from typing import DefaultDict, Counter as tCounter, Tuple
from .codfreq_types import NAPos, NAChar
from .posnas import iter_posnas

from .filename_helper import name_bamfile


GAP = ord(b'N')
ENCODING = 'UTF-8'


def sam2consensus(sampath, regions):
    nafreqs: DefaultDict[
        Tuple[NAPos, int],
        tCounter[NAChar]
    ] = defaultdict(Counter)
    for _, posnas in iter_posnas(sampath, progress=False):
        for refpos, i, na, _ in posnas:
            nafreqs[(refpos, i)][na] += 1
    nacons_lookup = {
        pos: nas.most_common(1)[0][0]
        for pos, nas in nafreqs.items()
    }

    results = []
    for region in regions:
        name = region['name']
        refpos_start = region['refStart']
        refpos_end = region['refEnd']
        consarr = bytearray()
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
        results.append({
            'name': name,
            'refStart': refpos_start,
            'refEnd': refpos_end,
            'consensus': consarr.decode(ENCODING)
        })
    return results


def create_untrans_region_consensus(seqname, profile):
    results = []
    for fragment in profile['fragmentConfig']:
        if 'fromFragment' in fragment:
            continue
        refname = fragment['fragmentName']
        samfile = name_bamfile(seqname, refname)
        regions = [
            region
            for region in profile['sequenceAssemblyConfig']
            if region.get('fromFragment') == refname
        ]
        results.extend(sam2consensus(samfile, regions))
    with open('{}.untrans.json'.format(seqname), 'w') as fp:
        json.dump(results, fp)
