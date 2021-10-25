from typing import Optional, TypedDict, Tuple, List, Dict


class Sequence(TypedDict):
    header: Optional[str]
    sequence: str


class FragmentConfig(TypedDict):
    fragmentName: Optional[str]
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


RefAAs = Dict[int, str]

RefFragment = Tuple[str, FragmentConfig, List[str]]
