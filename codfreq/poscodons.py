#! /usr/bin/env python

import os
import ray
from tqdm import tqdm
from statistics import mean
from more_itertools import chunked
from collections import defaultdict, Counter

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 30))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))
SITE_QUALITY_CUTOFF = int(os.environ.get('SITE_QUALITY_CUTOFF', 25))

ERR_OK = 0b00
ERR_TOO_SHORT = 0b01
ERR_LOW_QUAL = 0b10


def extract_codons(seq, qua, aligned_pairs):

    # pre-filter
    err = ERR_OK
    if not seq or len(seq) < LENGTH_CUTOFF:
        err |= ERR_TOO_SHORT
    if not seq or mean(qua) < OVERALL_QUALITY_CUTOFF:
        err |= ERR_LOW_QUAL
    if err:
        return None, err

    codons = {}
    codonqs = defaultdict(list)
    prev_refpos = 0
    prev_seqpos = 0
    min_cdpos = 0xffffffff
    max_cdpos = 0
    for seqpos, refpos in aligned_pairs:
        if refpos is None:
            # insertion
            refpos = prev_refpos - 1
        if seqpos is None:
            # deletion
            seqpos = prev_seqpos - 1

        if refpos > -1:
            codonpos = refpos // 3 + 1
            codonbp = refpos % 3
            min_cdpos = min(min_cdpos, codonpos)
            nas = []
            for n, q in zip(seq[prev_seqpos:seqpos + 1],
                            qua[prev_seqpos:seqpos + 1]):
                if q < SITE_QUALITY_CUTOFF:
                    n = '*'
                nas.append(n)
                codonqs[codonpos].append(q)
            if codonpos not in codons:
                codons[codonpos] = ['-'] * 3
            codons[codonpos][codonbp] = ''.join(nas)
            max_cdpos = max(max_cdpos, codonpos)
        prev_seqpos = seqpos + 1
        prev_refpos = refpos + 1
    # if prev_seqpos > -1:
    #     codonpos = refpos // 3 + 1
    #     codonbp = refpos % 3
    #     if codonpos not in codons:
    #         codons[codonpos] = ['-'] * 3
    #     codons[codonpos][codonbp] = seq[prev_seqpos]
    #     codonqs[codonpos].append(qua[prev_seqpos])
    #     max_cdpos = codonpos
    results = []
    for cdpos in range(min_cdpos, max_cdpos + 1):
        codon = ''.join(codons.get(cdpos, ['---']))
        qs = codonqs[cdpos]
        if qs:
            meanq = sum(qs) / len(qs)
        else:
            meanq = 0
        codon = codon if codon == '---' else codon.replace('-', '')
        # codon = codon.replace('-', '')
        if len(codon) - codon.count('*') < 3:
            continue
        codon = codon.replace('*', 'N')
        results.append((cdpos, codon, meanq))

    if results:
        lastpos, codon, q = results[-1]
        # remove insertion at the end of sequence read
        results[-1] = (lastpos, codon[0:3], q)
    return results, err


@ray.remote
def batch_extract_codons(input_data):
    return [extract_codons(*args) for args in input_data]


def iter_poscodons(all_paired_reads,
                   site_quality_cutoff=SITE_QUALITY_CUTOFF):

    futures = []
    all_headers = []
    all_paired_reads = list(all_paired_reads)
    pbar = tqdm(total=len(all_paired_reads))
    for partial in chunked(all_paired_reads, 20000):
        pbar.update(1)
        headers = []
        input_data = []
        for header, pair in partial:
            pbar.update(1)
            for read in pair:
                headers.append(header)
                input_data.append((
                    read.query_sequence,
                    read.query_qualities,
                    read.get_aligned_pairs(False)
                ))
        all_headers.append(headers)
        futures.append(batch_extract_codons.remote(input_data))
    pbar.close()
    for headers, all_results in zip(all_headers, ray.get(futures)):
        for header, (results, err) in zip(headers, all_results):
            poscodons = defaultdict(Counter)
            if err & ERR_TOO_SHORT or err & ERR_LOW_QUAL:
                continue
            for refpos, codon, q in results:
                if q < site_quality_cutoff:
                    continue
                poscodons[refpos][codon] = q
            if poscodons:
                yield header, [
                    # pos, codon, qual
                    (pos, *counter.most_common(1)[0])
                    for pos, counter in sorted(poscodons.items())
                ]