# cython: profile=True

from collections import Counter
from postalign.models.sequence import Sequence, PositionalSeqStr
from itertools import groupby

from .paired_reads import iter_paired_reads
from .posnas import iter_posnas

ENCODING = 'UTF-8'
GAP = ord(b'-')


def get_consensus_na(posnas):
    posnas = tuple(posnas)
    ((_, idx), cons_na), count = posnas[0]
    if idx == 0:
        return cons_na
    else:
        total = sum(count for _, count in posnas)
        if count > total * 0.9:
            return cons_na
        return None


def get_most_common_nas(sampath):
    nafreqs = Counter()
    all_paired_reads = iter_paired_reads(sampath)

    # flatten posnas generator (remove headers)
    nafreqs.update(
        ((pos, idx), na)
        for _, posnas in iter_posnas(all_paired_reads)
        for (pos, idx), na, _ in posnas
    )
    return {
        posidx: get_consensus_na(items)
        for posidx, items in groupby(
            sorted(
                nafreqs.items(),
                key=lambda item: (
                    item[0][0],  # posidx
                    -item[1]     # -count
                )
            ),
            lambda item: item[0][0]  # posidx
        )
        if posidx[1] == 0  # exclude all insertions when building consensus
    }


def iter_alignment_pairs(refseq_bytes, most_common_nas):
    for refpos0, refna in enumerate(refseq_bytes):
        refpos = refpos0 + 1
        if (refpos, 0) in most_common_nas:
            na = most_common_nas[(refpos, 0)]
            yield refna, na
            # i = 1
            # while (refpos, i) in most_common_nas:
            #     na = most_common_nas[(refpos, i)]
            #     if not na:
            #         break
            #     i += 1
            #     yield GAP, na
        else:
            yield refna, GAP


def assemble_alignment(refseq, most_common_nas):
    refseq_bytes = bytearray(refseq, ENCODING)
    return zip(*iter_alignment_pairs(refseq_bytes, most_common_nas))


def calc_pairwise_consensus(sampath, refseq):
    most_common_nas = get_most_common_nas(sampath)

    (alignedref,
     alignedseq) = assemble_alignment(refseq, most_common_nas)

    return (
        Sequence(
            header='ref',
            description='',
            seqtext=PositionalSeqStr.init_from_nastring(
                bytes(alignedref).decode(ENCODING)
            ),
            seqid=1,
            seqtype='NA',
            abs_seqstart=0,
            skip_invalid=True
        ),
        Sequence(
            header=sampath,
            description='',
            seqtext=PositionalSeqStr.init_from_nastring(
                bytes(alignedseq).decode(ENCODING)
            ),
            seqid=2,
            seqtype='NA',
            abs_seqstart=0,
            skip_invalid=True
        )
    )
