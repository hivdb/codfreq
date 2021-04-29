# cython: profile=True

from tqdm import tqdm

ENCODING = 'UTF-8'
GAP = ord(b'-')


def iter_single_read_posnas(seq, qua, aligned_pairs):
    seqchars = bytearray(seq, ENCODING)

    prev_refpos = 0
    prev_seqpos0 = 0
    idx = 0
    for seqpos0, refpos0 in aligned_pairs:

        if refpos0 is None:
            # insertion
            refpos = prev_refpos
            idx += 1
        else:
            refpos = refpos0 + 1
            idx = 0
            prev_refpos = refpos

        if seqpos0 is None:
            # deletion
            n, q = GAP, qua[prev_seqpos0]
        else:
            n, q = seqchars[seqpos0], qua[seqpos0]
            prev_seqpos0 = seqpos0

        if refpos == 0:
            continue

        yield (refpos, idx), n, q


def iter_posnas(all_paired_reads,
                site_quality_cutoff=0):

    all_paired_reads = list(all_paired_reads)
    for header, pair in tqdm(all_paired_reads):
        for read in pair:
            results = iter_single_read_posnas(
                read.query_sequence,
                read.query_qualities,
                read.get_aligned_pairs(False)
            )
            if site_quality_cutoff > 0:
                results = (
                    (posidx, na, q)
                    for posidx, na, q in results
                    if q >= site_quality_cutoff
                )
            yield header, results
