import cython  # type: ignore

from typing import (
    Tuple,
    List,
    Dict,
    Any,
    Optional,
    Union,
    Literal
)

from .codfreq_types import (
    AAPos,
    Header,
    Profile,
    NAPosRange,
    FragmentConfig,
    MainFragmentConfig,
    DerivedFragmentConfig,
    CodonAlignmentConfig
)
from .sam2codfreq_types import (
    TypedRefFragment
)

from .segfreq import DEFAULT_SEGMENT_SIZE, DEFAULT_SEGMENT_STEP


@cython.cfunc
@cython.inline
@cython.returns(list)
def get_ref_ranges(config: Dict) -> List[NAPosRange]:
    refstart = config.get('refStart')
    refend = config.get('refEnd')
    orig_refranges = config.get('refRanges')
    refranges: List[NAPosRange] = []
    if isinstance(orig_refranges, list):
        refranges = [(start, end) for start, end in orig_refranges]
    elif isinstance(refstart, int) and isinstance(refend, int):
        refranges = [(refstart, refend)]
    return refranges


@cython.ccall
@cython.inline
@cython.returns(list)
def get_ref_fragments(
    profile: Profile
) -> List[Tuple[
    Header,
    MainFragmentConfig,
    List[DerivedFragmentConfig]
]]:
    refname: str
    config: FragmentConfig
    ref_fragments: Dict[Header, TypedRefFragment] = {}
    frag_size_lookup: Dict[Header, AAPos] = {}
    for config in profile['fragmentConfig']:
        refname = config['fragmentName']
        refseq = config.get('refSequence')
        segment_size = config.get('segmentSize', DEFAULT_SEGMENT_SIZE)
        segment_step = config.get('segmentStep', DEFAULT_SEGMENT_STEP)
        if not isinstance(segment_size, int):
            raise TypeError('segmentSize must be an integer')
        if not isinstance(segment_step, int):
            raise TypeError('segmentStep must be an integer')
        if isinstance(refseq, str):
            ref_fragments[refname] = {
                'ref': {
                    'fragmentName': refname,
                    'refSequence': refseq,
                    'segmentSize': segment_size,
                    'segmentStep': segment_step
                },
                'fragments': []
            }

    for config in profile['fragmentConfig']:
        refname = config['fragmentName']
        fromref = config.get('fromFragment')
        gene = config.get('geneName')
        frequency_type = config.get('frequencyType', 'codon')
        refranges = get_ref_ranges(config)
        codon_alignment_raw: Any = config.get('codonAlignment')
        codon_alignment: Optional[Union[
            Literal[False], List[CodonAlignmentConfig]
        ]] = None
        if isinstance(codon_alignment_raw, list):
            codon_alignment = []
            for one in codon_alignment_raw:
                cda: CodonAlignmentConfig = {}
                if 'relRefStart' in one:
                    cda['relRefStart'] = one['relRefStart']
                if 'relRefEnd' in one:
                    cda['relRefEnd'] = one['relRefEnd']
                if 'windowSize' in one:
                    cda['windowSize'] = one['windowSize']
                if 'minGapDistance' in one:
                    cda['minGapDistance'] = one['minGapDistance']
                if 'relGapPlacementScore' in one:
                    cda['relGapPlacementScore'] = one['relGapPlacementScore']
                codon_alignment.append(cda)
        elif codon_alignment_raw is False:
            codon_alignment = False

        if frequency_type != 'codon' and \
                frequency_type != 'nucleotide':
            raise ValueError(
                'frequencyType must be "codon" or "nucleotide"')

        if (
            isinstance(fromref, str) and
            (gene is None or isinstance(gene, str)) and
            refranges
        ):

            ref_fragments[fromref]['fragments'].append({
                'fragmentName': refname,
                'fromFragment': fromref,
                'frequencyType': frequency_type,
                'geneName': gene,
                'refRanges': refranges,
                'codonAlignment': codon_alignment
            })
            frag_size_lookup[refname] = sum(
                (end - start + 1) for start, end in refranges
            ) // 3

    return [
        (refname, pair['ref'], pair['fragments'])
        for refname, pair in ref_fragments.items()
    ]
