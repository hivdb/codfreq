import typing
from typing import (
    DefaultDict,
    Tuple,
    List,
    TypedDict,
)
from .codfreq_types import (
    MainFragmentConfig,
    DerivedFragmentConfig,
)


class TypedRefFragment(TypedDict):
    ref: MainFragmentConfig
    genes: List[DerivedFragmentConfig]


CodonStat = DefaultDict[
    Tuple[str, int],
    typing.Counter[Tuple[bytes, ...]]
]

Qualities = typing.Counter[
    Tuple[str, int, Tuple[bytes, ...]]
]
