from typing import (
    TypedDict
)
from collections import Counter
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
    fragments: list[DerivedFragmentConfig]


CodonCounter = Counter[
    tuple[Header, AAPos, CodonText]
]

CodonCounterByFragPos = dict[
    tuple[Header, AAPos],
    Counter[CodonText]
]

FragmentGeneLookup = dict[
    Header, list[
        tuple[GeneText, AAPos]
        #                 ^
        #              AAOffset
    ]
]
