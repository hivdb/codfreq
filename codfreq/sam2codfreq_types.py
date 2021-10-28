from typing import (
    Tuple,
    List,
    TypedDict,
    DefaultDict,
    Counter
)
from .codfreq_types import (
    MainFragmentConfig,
    DerivedFragmentConfig,
)


class TypedRefFragment(TypedDict):
    ref: MainFragmentConfig
    genes: List[DerivedFragmentConfig]


CodonCounter = DefaultDict[
    Tuple[str, int],
    Counter[bytes]
]

QualityCounter = Counter[
    Tuple[str, int, bytes]
]
