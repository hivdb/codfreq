from typing import Optional, TypedDict, Tuple, List, Dict

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
    pair: Tuple[Optional[FASTQFileName], ...]
    n: int


class Sequence(TypedDict):
    header: Header
    sequence: SeqText


class MainFragmentConfig(TypedDict):
    fragmentName: Header
    refSequence: SeqText


class DerivedFragmentConfig(TypedDict):
    fragmentName: Header
    fromFragment: Header
    geneName: GeneText
    refStart: NAPos
    refEnd: NAPos


class FragmentConfig(TypedDict):
    fragmentName: str
    refSequence: Optional[str]
    fromFragment: Optional[str]
    geneName: Optional[str]
    refStart: Optional[int]
    refEnd: Optional[int]


class SequenceAssemblyConfig(TypedDict):
    name: Optional[str]
    geneName: Optional[str]
    fromFragment: Optional[str]
    refStart: Optional[int]
    refEnd: Optional[int]


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


RefAAs = Dict[AAPos, MultiAAText]
