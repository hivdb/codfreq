from typing import TypedDict, Literal
from pydantic import BaseModel

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


class MainFragmentConfig(BaseModel):
    """Reference fragment configuration.

    :param fragmentName: Name of the fragment.
    :type fragmentName: Header
    :param refSequence: Reference nucleotide sequence.
    :type refSequence: SeqText | None
    """

    fragmentName: Header
    refSequence: SeqText | None = None


class CodonAlignmentConfig(BaseModel):
    """Configuration options for codon alignment.

    :param relRefStart: Start position relative to the reference.
    :type relRefStart: NAPos | None
    :param relRefEnd: End position relative to the reference.
    :type relRefEnd: NAPos | None
    :param windowSize: Sliding window size in amino acids.
    :type windowSize: AAPos | None
    :param minGapDistance: Minimum nucleotide distance between gaps.
    :type minGapDistance: NAPos | None
    :param relGapPlacementScore: Relative gap placement score string.
    :type relGapPlacementScore: str | None
    """

    relRefStart: NAPos | None = None
    relRefEnd: NAPos | None = None
    windowSize: AAPos | None = None
    minGapDistance: NAPos | None = None
    relGapPlacementScore: str | None = None


class DerivedFragmentConfig(BaseModel):
    """Fragment derived from a reference fragment.

    :param fragmentName: Name of the fragment.
    :type fragmentName: Header
    :param fromFragment: Source reference fragment name.
    :type fromFragment: Header
    :param geneName: Gene identifier if any.
    :type geneName: GeneText | None
    :param refRanges: Reference coordinate ranges.
    :type refRanges: list[NAPosRange] | None
    :param codonAlignment: Codon alignment configuration or ``False`` to
        disable alignment.
    :type codonAlignment: None | (Literal[False] | list[CodonAlignmentConfig])
    """

    fragmentName: Header
    fromFragment: Header
    geneName: GeneText | None = None
    refRanges: list[NAPosRange] | None = None
    codonAlignment: None | (Literal[False] | list[CodonAlignmentConfig]) = None


FragmentConfig = MainFragmentConfig | DerivedFragmentConfig


class SequenceAssemblyConfig(BaseModel):
    """Configuration for assembling sequences from fragments.

    :param name: Name of the assembly region.
    :type name: str | None
    :param geneName: Associated gene name.
    :type geneName: str | None
    :param fromFragment: Source fragment name.
    :type fromFragment: str | None
    :param refStart: Start position in reference coordinates.
    :type refStart: int | None
    :param refEnd: End position in reference coordinates.
    :type refEnd: int | None
    """

    name: str | None = None
    geneName: str | None = None
    fromFragment: str | None = None
    refStart: int | None = None
    refEnd: int | None = None


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


class Profile(BaseModel):
    """Top level profile configuration.

    :param version: Profile schema version.
    :type version: str
    :param fragmentConfig: Fragment configuration list.
    :type fragmentConfig: list[FragmentConfig]
    :param sequenceAssemblyConfig: Assembly region configuration list.
    :type sequenceAssemblyConfig: list[SequenceAssemblyConfig]
    """

    version: str = ''
    fragmentConfig: list[FragmentConfig] = []
    sequenceAssemblyConfig: list[SequenceAssemblyConfig] = []


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
