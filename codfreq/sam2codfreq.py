from more_itertools import flatten
from collections import defaultdict, Counter

from .paired_reads import iter_paired_reads
from .poscodons import iter_poscodons
from .codonalign_consensus import codonalign_consensus
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


def sam2codfreq(
    samfile,
    ref,
    genes,
    site_quality_cutoff=0,
    log_format='text',
    include_partial_codons=False,
    **extras
):
    """Iterate CodFreq rows from a SAM/BAM file

    :param samfile: str of the SAM/BAM file path
    :param ref: dict of the reference configuration
    :param genes: list of dict of the gene configurations
    :param site_quality_cutoff: phred-score quality cutoff of each codon
                                position
    :param log_format: 'json' or 'text', default to 'text'
    :param include_partial_codons: if true output partial codons also; only for
                                   debug purpose
    :param **extras: any other variables to pass to the log method
    :return: CodFreq rows
    """

    codonstat = defaultdict(Counter)
    qualities = Counter()
    all_paired_reads = iter_paired_reads(samfile)

    # Iterate through the whole SAM/BAM and stat each individual gene position
    # and codon
    for _, poscodons in iter_poscodons(
        all_paired_reads,
        ref,
        genes,
        description=samfile,
        site_quality_cutoff=site_quality_cutoff,
        log_format=log_format,
        include_partial_codons=include_partial_codons,
        **extras
    ):
        for gene, refpos, codon, qua in poscodons:
            codonstat[(gene, refpos)][codon] += 1
            qualities[(gene, refpos, codon)] += qua
        del poscodons

    codonfreq = iter_codonfreq(codonstat, qualities)

    # Apply codon-alignment to consensus codons. This step is an imperfect
    # approach to address the out-frame deletions/insertions. However, it is
    # faster than codon-align each single read, since the process is right
    # now very slow and time-consuming.
    # This approach may need to be changed in the future when optimization was
    # done for the post-align codon-alignment program.
    yield from codonalign_consensus(codonfreq, ref, genes)


def sam2codfreq_all(
    name,
    fnpair,
    profile,
    site_quality_cutoff=0,
    log_format='text',
    include_partial_codons=False
):
    for refname, ref, genes in iter_ref_fragments(profile):
        samfile = name_samfile(name, refname)
        yield from sam2codfreq(
            samfile,
            ref,
            genes,
            site_quality_cutoff=site_quality_cutoff,
            log_format=log_format,
            include_partial_codons=include_partial_codons,
            fastqs=fnpair
        )
