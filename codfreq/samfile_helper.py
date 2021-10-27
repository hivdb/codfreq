import pysam  # type: ignore
import cython  # type: ignore
from typing import List, Tuple


@cython.ccall
@cython.inline
@cython.returns(list)
def chunked_samfile(
    samfile: str,
    chunk_size: int = 50000
) -> List[Tuple[int, int]]:
    """Find out positions in sam/bam file for given chunk_size"""

    chunks: List[Tuple[int, int]] = []
    cur_chunk_size: int
    cur_begin: int
    cur_end: int

    with pysam.AlignmentFile(samfile, 'rb') as samfp:
        cur_begin = samfp.tell()
        cur_chunk_size = 0

        # find out chunk start and end
        for _ in samfp:
            cur_chunk_size += 1
            if cur_chunk_size == chunk_size:
                cur_end = samfp.tell()
                chunks.append((cur_begin, cur_end))
                cur_begin = cur_end
                cur_chunk_size = 0
        cur_end = samfp.tell()
        if cur_end > cur_begin:
            chunks.append((cur_begin, cur_end))

    return chunks
