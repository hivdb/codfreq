import pysam  # type: ignore
import cython  # type: ignore
from tqdm import tqdm  # type: ignore
from collections import Counter
from more_itertools import unique_everseen

from typing import Any
from concurrent.futures import ProcessPoolExecutor

from .codfreq_types import (
    Header,
    AAPos,
    Profile,
    GeneText,
    CodFreqRow,
    CodonText,
    FragmentInterval,
    FASTQFileName,
    FragmentConfig,
    MainFragmentConfig,
    DerivedFragmentConfig
)
from .sam2codfreq_types import (
    CodonCounter,
    CodonCounterByFragPos,
    TypedRefFragment,
    FragmentGeneLookup
)
from .samfile_helper import chunked_samfile
from .json_progress import JsonProgress
from .poscodons import iter_poscodons, PosCodon
from .codonalign_consensus import codonalign_consensus
from .filename_helper import name_bamfile


CODFREQ_HEADER: list[str] = [
    'gene', 'position',
    'total', 'codon', 'count',
    'total_quality_score'
]


ENCODING: str = 'UTF-8'


@cython.cfunc
@cython.inline
@cython.returns(list)
def build_fragment_intervals(
    fragments: list[DerivedFragmentConfig]
) -> list[FragmentInterval]:
    """Extract fragment reference intervals for codon processing.

    :param fragments: Fragment configurations to convert.
    :type fragments: list[DerivedFragmentConfig]
    :returns: Intervals paired with fragment names.
    :rtype: list[FragmentInterval]
    """

    return [
        (fragment.refRanges, fragment.fragmentName)
        for fragment in fragments
    ]


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def get_ref_fragments(
    profile: Profile
) -> tuple[
    list[tuple[
        Header,
        MainFragmentConfig,
        list[DerivedFragmentConfig]
    ]],
    FragmentGeneLookup
]:
    refname: str
    config: FragmentConfig
    ref_fragments: dict[Header, TypedRefFragment] = {}
    frag_size_lookup: dict[Header, AAPos] = {}
    for config in profile.fragmentConfig:
        if isinstance(config, MainFragmentConfig):
            ref_fragments[config.fragmentName] = {
                'ref': config,
                'fragments': []
            }
    for config in profile.fragmentConfig:
        if not isinstance(config, DerivedFragmentConfig):
            continue
        refranges = config.refRanges
        ref_fragments[config.fromFragment]['fragments'].append(config)
        frag_size_lookup[config.fragmentName] = sum(
            (end - start + 1) for start, end in refranges
        ) // 3

    # build frag_gene_lookup
    gene_offsets: dict[GeneText, AAPos] = {}
    frag_gene_lookup: FragmentGeneLookup = {}
    for config in profile.fragmentConfig:
        gene = getattr(config, 'geneName', None)
        if not isinstance(gene, str):
            continue
        refname = config.fragmentName
        if gene not in gene_offsets:
            gene_offsets[gene] = 0
        if refname not in frag_gene_lookup:
            frag_gene_lookup[refname] = []
        frag_gene_lookup[refname].append((gene, gene_offsets[gene]))
        gene_offsets[gene] += frag_size_lookup[refname]

    return [
        (refname, pair['ref'], pair['fragments'])
        for refname, pair in ref_fragments.items()
    ], frag_gene_lookup


@cython.cfunc
@cython.inline
@cython.returns(dict)
def to_codon_counter_by_fragpos(
    codon_counter: CodonCounter
) -> CodonCounterByFragPos:
    """Reindex codon counts keyed by fragment and position."""

    fragment_name: Header
    refpos: AAPos
    codon: CodonText
    reformed: dict[
        tuple[Header, AAPos],
        Counter[CodonText]
    ] = {}

    for (fragment_name, refpos, codon), count in codon_counter.items():
        key = (fragment_name, refpos)
        if key not in reformed:
            reformed[key] = Counter()
        reformed[key][codon] = count

    return dict(reformed)


@cython.cfunc
@cython.inline
@cython.returns(list)
def get_codonfreq(
    codonstat_by_fragpos: CodonCounterByFragPos,
    qualities_by_fragpos: CodonCounterByFragPos,
    frag_gene_lookup: FragmentGeneLookup
) -> list[CodFreqRow]:
    """Calculate codon frequencies with quality scores."""

    fragment_name: Header
    gene: GeneText
    gene_offset: AAPos
    refpos: AAPos
    codon: CodonText
    count: int
    qua: int
    total: int
    rows: list[CodFreqRow] = []
    ordered_genes: list[GeneText] = list(unique_everseen([
        gene for genes in frag_gene_lookup.values() for gene, _ in genes
    ]))

    for (fragment_name, refpos), codons in codonstat_by_fragpos.items():
        total = sum(codons.values())
        for codon, count in codons.items():
            qua = qualities_by_fragpos[(fragment_name, refpos)][codon]
            for gene, gene_offset in frag_gene_lookup[fragment_name]:
                rows.append({
                    'gene': gene,
                    'position': refpos + gene_offset,
                    'total': total,
                    'codon': codon,
                    'count': count,
                    'total_quality_score': round(qua, 2)
                })
    rows.sort(key=lambda row: (
        ordered_genes.index(row['gene']),
        row['position'],
        row['codon']
    ))
    return rows


@cython.ccall
@cython.returns(tuple)
def sam2codfreq_between(
    samfile: str,
    samfile_start: int,
    samfile_end: int,
    fragment_intervals: list[FragmentInterval],
    site_quality_cutoff: int = 0
) -> tuple[CodonCounter, CodonCounter, int]:
    """subprocess function to call iter_poscodons and count codons

    The results of PosCodons are aggregated into codon counters to reduce
    memory usage before they are sent to the main process. This also reduces
    the computational complexity of the main process and relieves the risk of
    backpressure, which may lead to an out-of-memory issue.

    The out-of-memory issue was misidentified as a mysterious memory leak
    issue. It got worse when Cython was enabled. It turns out that Cython only
    speeded up the subprocesses and sent data to main process in shorter times.
    The main process was simply not fast enough to process such amount of huge
    data. The unprocessed data was piled in the memory of main process and
    caused out-of-memory eventually.

    The solution is simple, just processing the partial results in subprocess
    instead of in main process.
    """
    poscodons: list[PosCodon]
    fragment_name: str
    refpos: int
    codon: CodonText
    qua: int
    codonstat: CodonCounter = Counter()
    qualities: CodonCounter = Counter()
    num_row: int = 0

    for _, poscodons in iter_poscodons(
        samfile,
        samfile_start,
        samfile_end,
        fragment_intervals,
        site_quality_cutoff
    ):
        num_row += 1
        for fragment_name, refpos, codon, qua in poscodons:
            codonstat[(fragment_name, refpos, codon)] += 1
            qualities[(fragment_name, refpos, codon)] += qua

    return codonstat, qualities, num_row


def sam2codfreq(
    samfile: str,
    ref: MainFragmentConfig,
    fragments: list[DerivedFragmentConfig],
    workers: int,
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    include_partial_codons: bool = False,
    chunk_size: int = 25000,
    **extras: Any
) -> tuple[CodonCounterByFragPos, CodonCounterByFragPos]:
    """Return CodFreq rows from a SAM/BAM file.

    This function utilizes subprocesses to process alignment data from a
    segment of a SAM/BAM file and to aggregate the results into a codon counter
    and a quality score counter.  A progress bar is emitted via either ``tqdm``
    or :class:`~codfreq.json_progress.JsonProgress` depending on
    ``log_format`` and is updated as each chunk is processed.

    The codon counter and quality score counter are converted into CodFreq
    results. The consensus codon from CodFreq results are then processed by the
    codon-alignment function provided by ``hivdb/postalign``.

    :param samfile: str of the SAM/BAM file path
    :param ref: dict of the reference configuration
    :param fragments: list of dict of the fragment/gene configurations
    :param site_quality_cutoff: phred-score quality cutoff of each codon
                                position
    :param log_format: 'json' or 'text', default to 'text'
    :param include_partial_codons: if true output partial codons also; only for
                                   debug purpose
    :param **extras: any other variables to pass to the log method
    :return: CodFreq rows
    """

    total: int
    pbar: JsonProgress | tqdm | None

    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        total = samfp.mapped
        if log_format == 'json':
            pbar = JsonProgress(total=total, description=samfile, **extras)
        elif log_format == 'text':
            pbar = tqdm(total=total)
            pbar.set_description(f'Processing {samfile}')

    chunks: list[tuple[int, int]] = chunked_samfile(samfile, chunk_size)
    fragment_intervals: list[
        FragmentInterval
    ] = build_fragment_intervals(fragments)
    codonstat: CodonCounter = Counter()
    qualities: CodonCounter = Counter()

    with ProcessPoolExecutor(workers) as executor:

        for partial_codonstat, partial_qualities, num_row in executor.map(
            sam2codfreq_between,
            *zip(*[
                (
                    samfile,
                    samfile_begin,
                    samfile_end,
                    fragment_intervals,
                    site_quality_cutoff,
                )
                for samfile_begin, samfile_end in chunks
            ])
        ):
            codonstat.update(partial_codonstat)
            qualities.update(partial_qualities)
            if pbar:
                pbar.update(num_row)
        if pbar:
            pbar.close()

    codonstat_by_fragpos: CodonCounterByFragPos = (
        to_codon_counter_by_fragpos(codonstat)
    )
    qualities_by_fragpos: CodonCounterByFragPos = (
        to_codon_counter_by_fragpos(qualities)
    )

    # Apply codon-alignment to consensus codons. This step is an imperfect
    # approach to address the out-frame deletions/insertions. However, it
    # is faster than codon-align each single read, since the process is
    # right now very slow and time-consuming.  This approach may need to be
    # changed in the future when optimization was done for the post-align
    # codon-alignment program.
    codonstat_by_fragpos, qualities_by_fragpos = codonalign_consensus(
        codonstat_by_fragpos,
        qualities_by_fragpos,
        ref,
        fragments
    )

    return codonstat_by_fragpos, qualities_by_fragpos


def sam2codfreq_all(
    name: str,
    fnpair: tuple[FASTQFileName | None, ...],
    profile: Profile,
    workers: int,
    site_quality_cutoff: int = 0,
    log_format: str = 'text',
    include_partial_codons: bool = False
) -> list[CodFreqRow]:
    refname: str
    ref: MainFragmentConfig
    fragments: list[DerivedFragmentConfig]
    ref_fragments, frag_gene_lookup = get_ref_fragments(profile)
    all_codonstat_by_fragpos: CodonCounterByFragPos = {}
    all_qualities_by_fragpos: CodonCounterByFragPos = {}
    for refname, ref, fragments in ref_fragments:
        samfile: str = name_bamfile(name, refname, is_trimmed=True)
        codonstat_by_fragpos, qualities_by_fragpos = sam2codfreq(
            samfile,
            ref,
            fragments,
            workers=workers,
            site_quality_cutoff=site_quality_cutoff,
            log_format=log_format,
            include_partial_codons=include_partial_codons,
            fastqs=fnpair
        )
        all_codonstat_by_fragpos.update(codonstat_by_fragpos)
        all_qualities_by_fragpos.update(qualities_by_fragpos)

    codfreq_rows: list[CodFreqRow] = get_codonfreq(
        all_codonstat_by_fragpos,
        all_qualities_by_fragpos,
        frag_gene_lookup
    )
    return codfreq_rows
