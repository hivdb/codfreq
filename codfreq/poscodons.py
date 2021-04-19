#! /usr/bin/env python

import os
# import ray
from tqdm import tqdm
from statistics import mean
from collections import defaultdict, Counter

from .json_progress import JsonProgress

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 15))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))

ERR_OK = 0b00
ERR_TOO_SHORT = 0b01
ERR_LOW_QUAL = 0b10


def _extract_codons(
    read,
    ref_offset, site_quality_cutoff
):
    """
    Extract codons from samtools read object

    :param read: the samtools read object
    :param ref_offset: reference position offset specified by
        --reference-start
    :param site_quality_cutoff: site Phred score cutoff specified by
        --site-quality-cutoff
    """

    codons = {}
    codonqs = defaultdict(list)
    min_cdpos = 0xffffffff
    max_cdpos = 0

    seq = read.query_sequence
    qua = read.query_qualities
    aligned_pairs = read.get_aligned_pairs(False)

    # pre-filter
    err = ERR_OK
    len_seq = len(seq) if seq else -1
    mean_qua = mean(qua) if qua else 0
    if len_seq < LENGTH_CUTOFF:
        err |= ERR_TOO_SHORT
    if mean_qua < OVERALL_QUALITY_CUTOFF:
        err |= ERR_LOW_QUAL
    if err:
        return (codons, codonqs, min_cdpos, max_cdpos,
                err, len_seq, mean_qua)

    refpos_start = 0
    seqpos_start = 0
    for seqpos_end, refpos_end in aligned_pairs:
        if refpos_end:
            refpos_end -= ref_offset
        if refpos_end is None:
            # insertion
            refpos_end = refpos_start - 1
        if seqpos_end is None:
            # deletion
            seqpos_end = seqpos_start - 1

        if refpos_end > -1:
            codonpos = refpos_end // 3 + 1
            min_cdpos = min(min_cdpos, codonpos)
            nas = []
            for n, q in zip(seq[seqpos_start:seqpos_end + 1],
                            qua[seqpos_start:seqpos_end + 1]):
                if q < site_quality_cutoff:
                    n = '*'
                nas.append(n)
                codonqs[codonpos].append(q)
            if codonpos not in codons:
                codons[codonpos] = []
            if not nas:
                # deletion
                nas = ['-' * (refpos_end - refpos_start + 1)]
            codons[codonpos].append(''.join(nas))
            max_cdpos = max(max_cdpos, codonpos)
        seqpos_start = seqpos_end + 1
        refpos_start = refpos_end + 1
    return (codons, codonqs, min_cdpos, max_cdpos,
            err, len_seq, mean_qua)


def find_boundary_partial_codons(
    read,
    ref_offset, site_quality_cutoff
):
    """
    Extract codons from samtools read object

    :param read: the samtools read object
    :param ref_offset: reference position offset specified by
        --reference-start
    :param site_quality_cutoff: site Phred score cutoff specified by
        --site-quality-cutoff
    """
    (codons,
     codonqs,
     min_cdpos,
     max_cdpos,
     err,
     len_seq,
     mean_qua) = _extract_codons(read, ref_offset, site_quality_cutoff)

    results = []
    for cdpos in range(min_cdpos, max_cdpos + 1):
        codon = ''.join(codons.get(cdpos, ['---']))
        qs = codonqs[cdpos]
        if qs:
            meanq = sum(qs) / len(qs)
        else:
            meanq = 0
        if (
            (cdpos == min_cdpos and (
                len(codon) < 3 or codon[:1] in ('-', '*')
            )) or
            (cdpos == max_cdpos and (
                len(codon) < 3 or codon[-1:] in ('-', '*')
            ))
        ):
            # the left-most/right-most codon is a partial codon, e.g. --A/A--
            codon = codon.replace('*', 'N')
            results.append((cdpos, codon, meanq))
        else:
            continue
    return results, err, len_seq, mean_qua


def extract_codons(
    read,
    ref_offset, site_quality_cutoff,
    keep_outframe=True,
    include_boundary_partial_codons=False
):
    """
    Extract codons from samtools read object

    :param read: the samtools read object
    :param ref_offset: reference position offset specified by
        --reference-start
    :param site_quality_cutoff: site Phred score cutoff specified by
        --site-quality-cutoff
    :param keep_outframe: boolean option for keeping/dropping out-frame
        deletions; out-frame insertions are always preseved
    """
    (codons,
     codonqs,
     min_cdpos,
     max_cdpos,
     err,
     len_seq,
     mean_qua) = _extract_codons(read, ref_offset, site_quality_cutoff)

    results = []
    for cdpos in range(min_cdpos, max_cdpos + 1):
        is_boundary_partial_codon = False
        codon = ''.join(codons.get(cdpos, ['---']))
        qs = codonqs[cdpos]
        if qs:
            meanq = sum(qs) / len(qs)
        else:
            meanq = 0
        codon_wofs = codon if codon == '---' else codon.replace('-', '')
        if len(codon_wofs) < 3 and not keep_outframe:
            # drop out-frame deletions
            continue
        if cdpos == min_cdpos and (
            len(codon) < 3 or codon[:1] in ('-', '*')
        ):
            # the left-most codon is a partial codon, e.g. --A
            is_boundary_partial_codon = True
            if not include_boundary_partial_codons:
                continue
        elif cdpos == max_cdpos and (
            len(codon) < 3 or codon[-1:] in ('-', '*')
        ):
            # the right-most codon is a partial codon, e.g. A--
            is_boundary_partial_codon = True
            if not include_boundary_partial_codons:
                continue
        codon = codon.replace('*', 'N')
        results.append((cdpos, codon, meanq, is_boundary_partial_codon))

    if results:
        lastpos, codon, q, flag = results[-1]
        # remove insertion at the end of sequence read
        results[-1] = (lastpos, codon[0:3], q, flag)
    return results, err, len_seq, mean_qua


def iter_poscodons(all_paired_reads,
                   description,
                   reference_start=1,
                   fnpair=None,
                   site_quality_cutoff=0,
                   log_format='text',
                   include_boundary_partial_codons=False):

    all_paired_reads = list(all_paired_reads)
    total = len(all_paired_reads)
    if log_format == 'json':
        extras = {}
        if fnpair is not None:
            extras['fastqs'] = fnpair
        pbar = JsonProgress(
            total=total, description=description, **extras)
    else:
        pbar = tqdm(total=total)
        pbar.set_description('Processing {}'.format(description))

    for header, pair in all_paired_reads:
        pbar.update(1)
        for read in pair:
            results, err, len_seq, mean_qua = extract_codons(
                read,
                ref_offset=reference_start - 1,
                site_quality_cutoff=site_quality_cutoff,
                include_boundary_partial_codons=include_boundary_partial_codons
            )
            poscodons = defaultdict(Counter)
            if err & ERR_TOO_SHORT or err & ERR_LOW_QUAL:
                continue
            for refpos, codon, q, is_boundary_partial_codon in results:
                if q < site_quality_cutoff:
                    continue
                poscodons[refpos][(codon, is_boundary_partial_codon)] = q
            if poscodons:
                result = []
                for pos, counter in sorted(poscodons.items()):
                    ((codon, is_boundary_partial_codon),
                     qual), = counter.most_common(1)
                    result.append((pos, codon, qual,
                                   is_boundary_partial_codon))
                yield header, result
    pbar.close()
