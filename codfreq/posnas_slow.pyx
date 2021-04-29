# cython: profile=True
from cpython cimport array
import array

from tqdm import tqdm
from itertools import groupby


def iter_nqpairs(unsigned char[:] seq, unsigned char[:] qua, list aligned_pairs):

    cdef int prev_refpos = 0
    cdef int prev_seqpos0 = 0
    cdef int refpos
    cdef unsigned char n
    cdef unsigned char q
    for seqpos0, refpos0 in aligned_pairs:

        if refpos0 is None:
            # insertion
            refpos = prev_refpos
        else:
            refpos = refpos0 + 1
            prev_refpos = refpos

        if seqpos0 is None:
            # deletion
            n, q = b'-', qua[prev_seqpos0]
        else:
            n, q = seq[seqpos0], qua[seqpos0]
            prev_seqpos0 = seqpos0

        if refpos == 0:
            continue

        yield (refpos, n, q)


def get_nqpair_key(tuple nqpair):
    return nqpair[0]


def iter_single_read_posnas(unsigned char[:] seq, unsigned char[:] qua, list aligned_pairs):
    nqpairs = iter_nqpairs(seq, qua, aligned_pairs)
    for pos, nqs in groupby(nqpairs, get_nqpair_key):
        nqs = list(nqs)
        nqs_size = len(nqs)
        if nqs_size == 1:
            nas = array.array('B', [nqs[0][1]])
            meanq = nqs[0][2]
        else:
            naslist = [nq[1] for nq in nqs]
            nas = array.array('B', naslist)
            meanq = sum(nq[2] for nq in nqs) / nqs_size
        yield pos, nas, meanq


def iter_posnas(all_paired_reads,
                site_quality_cutoff=0):

    all_paired_reads = list(all_paired_reads)
    for header, pair in tqdm(all_paired_reads):
        for read in pair:
            results = iter_single_read_posnas(
                array.array('B', list(read.query_sequence.encode('UTF-8'))),
                read.query_qualities,
                read.get_aligned_pairs(False)
            )
            if site_quality_cutoff > 0:
                results = (
                    (pos, nas, q)
                    for pos, nas, q in results
                    if q >= site_quality_cutoff
                )
            yield header, results
