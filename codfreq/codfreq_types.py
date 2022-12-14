from typing import Optional, TypedDict, Tuple, List, Dict, Union, Literal

FASTQFileName = str
Header = str
SeqText = str
AAPos = int
NAPos = int
GeneText = str
NAChar = int
MultiNAText = bytes
CodonText = MultiNAText
AAChar = int
MultiAAText = bytes

NAPosRange = Tuple[NAPos, NAPos]


class PairedFASTQ(TypedDict):
    name: Header
    pair: Tuple[FASTQFileName, Optional[FASTQFileName]]
    n: int


class Sequence(TypedDict):
    header: Header
    sequence: SeqText


class MainFragmentConfig(TypedDict):
    fragmentName: Header
    refSequence: SeqText


class CodonAlignmentConfig(TypedDict, total=False):
    relRefStart: NAPos
    relRefEnd: NAPos
    windowSize: Optional[AAPos]
    minGapDistance: Optional[NAPos]
    relGapPlacementScore: Optional[str]


class DerivedFragmentConfig(TypedDict, total=False):
    fragmentName: Header
    fromFragment: Header
    refSequence: None
    geneName: Optional[GeneText]
    refRanges: List[NAPosRange]
    codonAlignment: Optional[Union[
        Literal[False], List[CodonAlignmentConfig]
    ]]


FragmentConfig = Union[
    MainFragmentConfig,
    DerivedFragmentConfig
]


class SequenceAssemblyConfig(TypedDict):
    name: Optional[str]
    geneName: Optional[str]
    fromFragment: Optional[str]
    refStart: Optional[int]
    refEnd: Optional[int]


class NARegionConfig(TypedDict):
    name: str
    fromFragment: str
    refStart: int
    refEnd: int


class RegionalConsensus(TypedDict):
    name: str
    refStart: NAPos
    refEnd: NAPos
    consensus: str


class Profile(TypedDict):
    fragmentConfig: List[FragmentConfig]
    sequenceAssemblyConfig: List[SequenceAssemblyConfig]


class CodFreqRow(TypedDict):
    gene: GeneText
    position: AAPos
    total: int
    codon: CodonText
    count: int
    total_quality_score: int


#                                 refStart refEnd
#                                     v      v
FragmentInterval = Tuple[List[Tuple[NAPos, NAPos]], Header]

RefAAs = Dict[AAPos, MultiAAText]
