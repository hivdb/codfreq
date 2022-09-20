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
    refStart: NAPos
    refEnd: NAPos
    windowSize: Optional[AAPos]
    minGapDistance: Optional[NAPos]
    gapPlacementScore: Optional[str]


class DerivedFragmentConfig(TypedDict, total=False):
    fragmentName: Header
    fromFragment: Header
    refSequence: None
    geneName: Optional[GeneText]
    refStart: NAPos
    refEnd: NAPos
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


#                      refStart refEnd
#                          v      v
FragmentInterval = Tuple[NAPos, NAPos, Header]

RefAAs = Dict[AAPos, MultiAAText]
