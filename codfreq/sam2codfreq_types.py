from typing import (
    Dict,
    Tuple,
    List,
    TypedDict,
    Counter
)
from .codfreq_types import (
    Header,
    AAPos,
    GeneText,
    CodonText,
    MainFragmentConfig,
    DerivedFragmentConfig
)


class TypedRefFragment(TypedDict):
    ref: MainFragmentConfig
    fragments: List[DerivedFragmentConfig]


CodonCounter = Counter[
    Tuple[Header, AAPos, CodonText]
]

CodonCounterByFragPos = Dict[
    Tuple[Header, AAPos],
    Counter[CodonText]
]

FragmentGeneLookup = Dict[
    Header, List[
        Tuple[GeneText, AAPos]
        #                 ^
        #              AAOffset
    ]
]
