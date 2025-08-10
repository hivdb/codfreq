from typing import TypedDict, Literal

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

NAPosRange = tuple[NAPos, NAPos]


class PairedFASTQ(TypedDict):
    name: Header
    pair: tuple[FASTQFileName, FASTQFileName | None]
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
    windowSize: AAPos | None
    minGapDistance: NAPos | None
    relGapPlacementScore: str | None


class DerivedFragmentConfig(TypedDict, total=False):
    fragmentName: Header
    fromFragment: Header
    refSequence: None
    geneName: GeneText | None
    refRanges: list[NAPosRange]
    codonAlignment: None | (
        Literal[False] | list[CodonAlignmentConfig]
    )


FragmentConfig = MainFragmentConfig | DerivedFragmentConfig


class SequenceAssemblyConfig(TypedDict):
    name: str | None
    geneName: str | None
    fromFragment: str | None
    refStart: int | None
    refEnd: int | None


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
    version: str
    fragmentConfig: list[FragmentConfig]
    sequenceAssemblyConfig: list[SequenceAssemblyConfig]


class CodFreqRow(TypedDict):
    gene: GeneText
    position: AAPos
    total: int
    codon: CodonText
    count: int
    total_quality_score: int


#                                 refStart refEnd
#                                     v      v
FragmentInterval = tuple[list[tuple[NAPos, NAPos]], Header]

RefAAs = dict[AAPos, MultiAAText]
