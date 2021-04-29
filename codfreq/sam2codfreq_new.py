from more_itertools import flatten
from collections import defaultdict, Counter

from .paired_reads import iter_paired_reads
from .poscodons_new import iter_poscodons
from .codonalign_consensus import codonalign_consensus
# from .codonalign_caller import codonalign_caller
from .filename_helper import (
    name_samfile
)

CODFREQ_HEADER = [
    'gene', 'position',
    'total', 'codon', 'count',
    'total_quality_score'
]

ENCODING = 'UTF-8'


def iter_ref_fragments(profile):
    ref_fragments = {}
    for config in profile['fragmentConfig']:
        if 'refSequence' not in config:
            continue
        refname = config['fragmentName']
        ref_fragments[refname] = {
            'ref': config,
            'genes': []
        }
    for config in profile['fragmentConfig']:
        if 'geneName' not in config:
            continue
        refname = config['fragmentName']
        if 'fromFragment' in config:
            refname = config['fromFragment']
        ref_fragments[refname]['genes'].append(config)
    for refname, pair in ref_fragments.items():
        yield refname, pair['ref'], pair['genes']


def build_remap_lookup(pre_refseq, pre_queryseq, post_refseq, post_queryseq):
    cur_idx = 0
    cur_refpos = None
    post_map = {}
    for refpos, querypos in zip(post_refseq.seqtext._seq_pos,
                                post_queryseq.seqtext._seq_pos):
        if refpos > 0:
            cur_refpos = refpos
            cur_idx = 0
        else:
            cur_idx += 1
        if querypos > 0:
            post_map[querypos] = (cur_refpos, cur_idx)

    cur_idx = 0
    cur_refpos = None
    lookup = {}
    for refpos, querypos in zip(pre_refseq.seqtext._seq_pos,
                                pre_queryseq.seqtext._seq_pos):
        if refpos > 0:
            cur_refpos = refpos
            cur_idx = 0
        else:
            cur_idx += 1
        if querypos > 0:
            preval = (cur_refpos, cur_idx)
            postval = post_map[querypos]
            if preval != postval:
                lookup[preval] = postval

    return lookup


def iter_codonfreq(codonstat, qualities):
    for (gene, refpos), codons in codonstat.items():
        total = sum(codons.values())
        for codon, count in sorted(codons.most_common()):
            qua = qualities[(gene, refpos, codon)]
            yield {
                'gene': gene,
                'position': refpos,
                'total': total,
                'codon': bytes(flatten(codon)).decode(ENCODING),
                'count': count,
                'total_quality_score': round(qua, 2)
            }


def sam2codfreq_single_ref(
    fnpair,
    pattern,
    ref,
    genes,
    site_quality_cutoff=0,
    log_format='text',
    include_partial_codons=False
):
    refname = ref['fragmentName']
    # refseq = ref['refSequence']
    samfile = name_samfile(fnpair, pattern, refname)
    # refseq, queryseq = calc_pairwise_consensus(samfile, refseq)
    # post_refseq, post_queryseq = codonalign_caller(refseq, queryseq, ref)
    # remap_lookup = build_remap_lookup(
    #     refseq, queryseq,
    #     post_refseq, post_queryseq)
    # print(lookup)
    # print(list(zip(refseq.seqtext._seq_pos, queryseq.seqtext._seq_pos)))
    # print(list(zip(post_refseq.seqtext._seq_pos,
    #                post_queryseq.seqtext._seq_pos)))
    # return

    codonstat = defaultdict(Counter)
    qualities = Counter()
    all_paired_reads = iter_paired_reads(samfile)

    for _, poscodons in iter_poscodons(
        all_paired_reads,
        {},  # remap_lookup,
        ref,
        genes,
        fnpair=fnpair,
        description=samfile,
        site_quality_cutoff=site_quality_cutoff,
        log_format=log_format,
        include_partial_codons=include_partial_codons
    ):
        for gene, refpos, codon, qua in poscodons:
            codonstat[(gene, refpos)][codon] += 1
            qualities[(gene, refpos, codon)] += qua
        del poscodons

    codonfreq = iter_codonfreq(codonstat, qualities)
    # yield from codonfreq
    yield from codonalign_consensus(codonfreq, ref, genes)


def sam2codfreq_new(
    fnpair,
    pattern,
    profile,
    site_quality_cutoff=0,
    log_format='text',
    include_partial_codons=False
):
    for refname, ref, genes in iter_ref_fragments(profile):
        yield from sam2codfreq_single_ref(
            fnpair,
            pattern,
            ref,
            genes,
            site_quality_cutoff=site_quality_cutoff,
            log_format=log_format,
            include_partial_codons=include_partial_codons
        )
