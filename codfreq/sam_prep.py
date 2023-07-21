from typing import Tuple, List, Counter as tCounter, Optional
from collections import Counter

from pysam import AlignmentFile


def squash_gaps(
    cigartuples: Tuple[Tuple[int, int], ...],
    max_squashing_distance: int = 10
) -> Tuple[Tuple[int, int], ...]:
    """Squash gaps in cigar tuples.

    This function squashs neighboring deletions and insertions (the distance of
    deletions and insertions are not greater than `max_squashing_distance`)
    into one operation.
    """
    # Squash neighboring deletions and insertions
    tmp_cigartuples = list(cigartuples)
    prev_indel_idx: Optional[int] = None
    prev_indel_dist: int = 0
    for idx, (op, length) in enumerate(tmp_cigartuples):
        if op in (1, 2):
            # op is I/D
            if prev_indel_idx is not None and \
                    prev_indel_dist < max_squashing_distance:
                # merge the current indel into the previous one
                prev_op, prev_length = tmp_cigartuples[prev_indel_idx]
                if op == prev_op:
                    prev_length += length
                else:
                    prev_length -= length
                if prev_length < 0:
                    prev_op = op
                    prev_length = -prev_length
                tmp_cigartuples[prev_indel_idx] = (prev_op, prev_length)
                tmp_cigartuples[idx] = (op, 0)
                # reset prev_indel_dist
                prev_indel_dist = 0
            else:
                prev_indel_idx = idx
                prev_indel_dist = 0
        else:
            prev_indel_dist += length

    # Remove zero-length operations and merge adjacent operations
    result_cigartuples: List[Tuple[int, int]] = []
    for op, length in tmp_cigartuples:
        if length == 0:
            continue
        if result_cigartuples and result_cigartuples[-1][0] == op:
            result_cigartuples[-1] = (op, result_cigartuples[-1][1] + length)
        else:
            result_cigartuples.append((op, length))
    return tuple(result_cigartuples)


def count_indel_positions(
    cigartuples: Tuple[Tuple[int, int], ...],
    ref_start: int,
    indel_counter: tCounter[int]
) -> None:
    """Count indel positions in cigar tuples.

    This function counts the positions where insertions and deletions occur.
    """
    cum_pos = ref_start
    for op, length in cigartuples:
        if op in (1, 2):
            # op is I/D
            indel_counter[cum_pos] += 1
        if op in (0, 2, 3, 7, 8):
            # op is M/D/N/=/X, refseq is consumed
            cum_pos += length


def prepare_sam(sam_file: str, out_file: str) -> None:
    """Prepare SAM/BAM alignment files for downstream analysis."""
    cigars: List[Tuple[Tuple[int, int], ...]] = []
    indel_counter: tCounter[int] = Counter()

    with AlignmentFile(sam_file, "r") as sam, \
            AlignmentFile(out_file, "w", template=sam) as out:
        fp = sam.tell()
        for read in sam:
            if read.is_unmapped:
                cigars.append(())
                continue
            cigars.append(squash_gaps(read.cigartuples))
            count_indel_positions(
                cigars[-1], read.reference_start, indel_counter)

        sam.seek(fp)
        for read, cigar in zip(sam, cigars):
            # TODO: modify cigar according to indel_counter
            read.cigartuples = cigar
            out.write(read)
