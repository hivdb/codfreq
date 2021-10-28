from typing import (
    Tuple,
    List,
    TypedDict,
    Counter
)
from .codfreq_types import (
    MainFragmentConfig,
    DerivedFragmentConfig,
)


class TypedRefFragment(TypedDict):
    ref: MainFragmentConfig
    genes: List[DerivedFragmentConfig]


CodonCounter = Counter[
    Tuple[str, int, bytes]
]
