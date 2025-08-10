from .codfreq_types import NAPos, NAChar


PosNA = tuple[
    NAPos,   # refpos
    int,     # insertion_index
    NAChar,  # na
    int      # qua
]
