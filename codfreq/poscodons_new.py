# cython: profile=True

import os
from tqdm import tqdm
from itertools import groupby
from operator import itemgetter

from .json_progress import JsonProgress
from .posnas import iter_single_read_posnas

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 15))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))

# frameshift types
REPEAT = 0x01
SKIP = 0x02

UNSEQ = b'.'


def build_gene_intervals(genes):
    return [(
        genedef['refStart'],
        genedef['refEnd'],
        genedef['geneName']
    ) for genedef in genes]


def group_posnas(posnas, gene_intervals):
    """Group posnas by position and add gene AA position"""
    grouped_posnas = groupby(posnas, key=lambda posna: posna[0][0])
    for gene_refstart, gene_refend, gene in gene_intervals:
        pos = 0
        while pos < gene_refend:
            try:
                pos, na_and_ins = next(grouped_posnas)
                if pos < gene_refstart:
                    continue
                aapos = (pos - gene_refstart) // 3 + 1
                na_and_ins_tuple = tuple(na_and_ins)
                # if gene == 'nsp1' and aapos == 90:
                #     print(gene, aapos, na_and_ins_tuple)

                yield (
                    gene,
                    aapos,
                    na_and_ins_tuple
                )
            except StopIteration:
                break


def find_overlapped_genes(gene_intervals, read_refstart, read_refend):
    for gene_refstart, gene_refend, gene in gene_intervals:
        if read_refend < gene_refstart:
            break
        if read_refstart > gene_refend:
            continue
        yield gene_refstart, gene_refend, gene


def get_comparable_codon(codon):
    codon_nas = tuple(
        tuple(na for _, na, _ in posnas)
        for _, _, posnas in codon
    )
    is_partial = len(codon_nas) < 3
    return codon_nas, is_partial


def posnas2poscodons(
    posnas,
    gene_intervals,
    read_refstart,  # 1-based first aligned refpos
    read_refend,    # 1-based last aligned refpos
    site_quality_cutoff
):
    genes = find_overlapped_genes(gene_intervals, read_refstart, read_refend)
    grouped_posnas = group_posnas(posnas, genes)
    for (gene, aapos), codon in groupby(grouped_posnas, key=itemgetter(0, 1)):
        codon = tuple(codon)
        meanq = [
            qua
            for _, _, posnas in codon
            for _, _, qua in posnas
        ]
        meanq = sum(meanq) / len(meanq) if meanq else 0
        if meanq < site_quality_cutoff:
            continue
        nas, is_partial = get_comparable_codon(codon)
        if is_partial:
            continue
        yield gene, aapos, nas, meanq


def iter_poscodons(
    all_paired_reads,
    remap_lookup,
    ref,
    genes,
    description,
    site_quality_cutoff=0,
    log_format='text',
    **extras
):
    gene_intervals = build_gene_intervals(genes)

    all_paired_reads = list(all_paired_reads)
    total = len(all_paired_reads)
    if log_format == 'json':
        pbar = JsonProgress(
            total=total, description=description, **extras)
    else:
        pbar = tqdm(total=total)
        pbar.set_description('Processing {}'.format(description))

    for header, pair in all_paired_reads:
        pbar.update(1)
        for read in pair:

            posnas = iter_single_read_posnas(
                read.query_sequence,
                read.query_qualities,
                read.get_aligned_pairs(False)
            )
            remapped_posnas = (
                (remap_lookup.get(orig_posidx, orig_posidx), n, q)
                for orig_posidx, n, q in posnas
            )

            poscodons = posnas2poscodons(
                remapped_posnas,
                gene_intervals,
                read.reference_start + 1,  # pysam has 0-based numbering
                read.reference_end,  # "reference_end points to one past the
                                     #  last aligned residue."
                site_quality_cutoff
            )
            yield header, poscodons
    pbar.close()
