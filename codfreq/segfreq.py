import cython  # type: ignore
import re
import csv
from collections import deque, Counter

from typing import (
    Sequence, Dict, Tuple, TextIO, List,
    Optional, Deque, Counter as tCounter
)

from more_itertools import pairwise

from .posnas import PosNA, merge_posnas

DEFAULT_SEGMENT_SIZE: int = 3
DEFAULT_SEGMENT_STEP: int = 1
DEFAULT_TOP_N_SEEDS: int = 10
DEFAULT_CONSENSUS_LEVEL: float = 1.


@cython.ccall
@cython.inline
@cython.returns(cython.long)
def get_segment_pos(segment: Tuple[Optional[PosNA], ...]) -> int:
    for idx, node in enumerate(segment):
        if node is not None:
            return node.pos - idx
    raise ValueError("Segment is malformed")


@cython.ccall
@cython.inline
@cython.returns(tuple)
def remove_first_n_pos(
    segment: Tuple[Optional[PosNA], ...],
    n: int
) -> Tuple[Optional[PosNA], ...]:
    n_count: int = 0
    for idx, node in enumerate(segment):
        if node is None or node.bp == 0:
            n_count += 1
        if n_count > n:
            return segment[idx:]
    raise ValueError("Segment is malformed")


@cython.ccall
@cython.inline
@cython.returns(tuple)
def remove_last_n_pos(
    segment: Tuple[Optional[PosNA], ...],
    n: int
) -> Tuple[Optional[PosNA], ...]:
    n_count: int = 0
    for idx, node in enumerate(reversed(segment)):
        if node is None or node.bp == 0:
            n_count += 1
        if n_count == n:
            return segment[:len(segment) - idx - 1]
    raise ValueError("Segment is malformed")


@cython.ccall
@cython.inline
@cython.returns(cython.bint)
def is_continuous(
    left_segment: Tuple[Optional[PosNA], ...],
    right_segment: Tuple[Optional[PosNA], ...],
    segment_step: int
) -> bool:
    result: bool = (
        remove_first_n_pos(left_segment, segment_step) ==
        remove_last_n_pos(right_segment, segment_step)
    )
    return result


@cython.ccall
@cython.inline
@cython.returns(tuple)
def mask_segment(
    segment: Tuple[Optional[PosNA], ...],
    min_pos: int,
    max_pos: int
) -> Tuple[Optional[PosNA], ...]:
    return tuple([
        node
        if node is None or min_pos <= node.pos <= max_pos else
        None
        for node in segment
    ])


@cython.cclass
class SegFreq:
    """SegFreq class to represent an NGS alignment."""

    segment_size: int = cython.declare(cython.int, visibility="public")
    segment_step: int = cython.declare(cython.int, visibility="public")
    _segments: Dict[
        int,
        tCounter[Tuple[Optional[PosNA], ...]]
    ] = cython.declare(dict, visibility="private")
    _max_segpos: int = cython.declare(cython.int, visibility='private')

    def __init__(self: 'SegFreq', segment_size: int, segment_step: int):
        if segment_step < 1:
            raise ValueError("Segment step must be at least 1")
        if segment_size - segment_step < 2:
            raise ValueError(
                "Segment size must be at least segment step + 2"
            )
        self.segment_size = segment_size
        self.segment_step = segment_step
        self._segments = {}
        self._max_segpos = 0

    @cython.ccall
    @cython.locals(
        na_size=cython.long
    )
    @cython.returns(dict)
    def get_frequency(
        self: 'SegFreq',
        positions: Sequence[int],
        na_size: int = 3
    ) -> Dict[bytes, int]:
        """Get the NAs combination counts for given positions."""
        pos: int = cython.declare(cython.long)
        accessed: bool = cython.declare(cython.bint)
        nas: bytearray = cython.declare(bytearray)
        nas_bytes: bytes = cython.declare(bytes)
        positions = list(positions)
        while len(positions) < na_size:
            positions.append(positions[-1] + 1)
        segpos: int = min(positions)
        segpos = segpos - (segpos - 1) % self.segment_step
        if segpos > self._max_segpos:
            segpos = self._max_segpos
        for pos in positions:
            if pos >= self._max_segpos + self.segment_size:
                # current segment doesn't contain the position, skip
                continue
            if pos < segpos or pos >= segpos + self.segment_size:
                raise ValueError('Positions are too far apart')
        nas_counts: tCounter[bytes] = Counter()
        for segment, count in self._segments.get(segpos, {}).items():
            nas = bytearray()
            for pos in positions:
                accessed = False
                for node in segment:
                    if node and node.pos == pos:
                        nas.append(node.na)
                        accessed = True
                    elif node and node.pos > pos:
                        break
                if not accessed:
                    break
            else:
                nas_bytes = bytes(nas)
                nas_counts[nas_bytes] += count
        return dict(nas_counts)

    @cython.ccall
    @cython.locals(
        pos=cython.long
    )
    @cython.returns(dict)
    def get_pos_nas(
        self: 'SegFreq',
        pos: int
    ) -> Dict[bytes, int]:
        """Get the nucleotide counts for a given position."""
        segpos: int = pos - (pos - 1) % self.segment_step
        if segpos > self._max_segpos:
            # handle the case where the position is beyond the last segment
            segpos = self._max_segpos
        na_counts: tCounter[bytes] = Counter()
        nas: bytearray = cython.declare(bytearray)
        for segment, count in self._segments.get(segpos, {}).items():
            nas = bytearray()
            for node in segment:
                if node and node.pos == pos:
                    nas.append(node.na)
            if nas:
                na_counts[bytes(nas)] += count
        return dict(na_counts)

    @cython.ccall
    @cython.locals(
        pos_start=cython.long,
        pos_end=cython.long
    )
    @cython.returns(tuple)
    def get_consensus(
        self: 'SegFreq',
        pos_start: int,
        pos_end: int,
        level: float = DEFAULT_CONSENSUS_LEVEL
    ) -> Tuple[Optional[PosNA], ...]:
        """Get consensus segment for a range of positions."""
        pos: int = cython.declare(cython.long)
        consensus: Dict[Tuple[int, int], tCounter[PosNA]] = {}
        real_pos_start: int = pos_start - (pos_start - 1) % self.segment_step
        real_pos_end: int = pos_end - (pos_end - 1) % self.segment_step
        pos_total: Dict[int, int] = {}

        for pos in range(real_pos_start, real_pos_end + 1, self.segment_step):
            pos_until: int = pos + self.segment_step
            if pos in self._segments:
                pos_segments = self._segments[pos]
                total = sum(pos_segments.values())
                for segment, count in pos_segments.most_common():
                    for node in segment:
                        if node is None or \
                                node.pos >= pos_until or \
                                node.pos < pos_start or node.pos > pos_end:
                            continue
                        pos_total[node.pos] = total
                        nodekey = (node.pos, node.bp)
                        if nodekey not in consensus:
                            consensus[nodekey] = Counter()
                        consensus[nodekey][node] += count
                    if level >= 1.:
                        # only record the most common NA's for 100%
                        break
        result: List[Optional[PosNA]] = []
        for pos in range(pos_start, pos_end + 1):
            bp = 0
            while (pos, bp) in consensus:
                posnas: tCounter[PosNA] = consensus.get((pos, bp), Counter())
                min_count: float = pos_total[pos] * level
                qualified: List[PosNA] = [
                    posna
                    for posna, count in posnas.most_common()
                    if level >= 1. or count >= min_count
                ]
                if qualified:
                    result.append(merge_posnas(*qualified))
                bp += 1
        return tuple(result)

    @cython.ccall
    @cython.locals(
        pos_start=cython.long,
        pos_end=cython.long,
        top_n_seeds=cython.int
    )
    @cython.returns(dict)
    def get_patterns(
        self: 'SegFreq',
        pos_start: int,
        pos_end: int,
        top_n_seeds: int = DEFAULT_TOP_N_SEEDS
    ) -> Dict[Tuple[PosNA, ...], Tuple[int, float]]:
        """Get segments between two positions."""
        real_pos_start = pos_start - (pos_start - 1) % self.segment_step
        real_pos_end = pos_end - self.segment_size + 1
        real_pos_end = real_pos_end + (1 - real_pos_end) % self.segment_step
        real_pos_end = (
            real_pos_start
            if real_pos_start > real_pos_end else
            real_pos_end
        )
        segments_between_pcnt: Dict[
            int,
            tCounter[Tuple[Optional[PosNA], ...]]
        ] = {}
        segments_between_ct: Dict[
            int,
            tCounter[Tuple[Optional[PosNA], ...]]
        ] = {}
        seed_segments_pcnt: tCounter[
            Tuple[int, Tuple[Optional[PosNA], ...]]
        ] = Counter()
        seed_segments_ct: tCounter[
            Tuple[int, Tuple[Optional[PosNA], ...]]
        ] = Counter()
        for pos in range(real_pos_start, 1 + real_pos_end, self.segment_step):
            segments_between_pcnt[pos] = Counter()
            segments_between_ct[pos] = Counter()
            total = sum(self._segments.get(pos, {}).values())
            for segment, count in self._segments.get(pos, {}).items():
                segment = mask_segment(segment, pos_start, pos_end)
                segments_between_pcnt[pos][segment] += count * 10000 // total
                segments_between_ct[pos][segment] += count
                seed_segments_pcnt[(pos, segment)] += count * 10000 // total
                seed_segments_ct[(pos, segment)] += count

        patterns_pcnt: tCounter[Tuple[PosNA, ...]] = Counter()
        patterns_ct: tCounter[Tuple[PosNA, ...]] = Counter()
        pattern_pcnt: int
        pattern_count: int
        seed_segment: Tuple[Optional[PosNA], ...]
        prev_segment: Tuple[Optional[PosNA], ...]

        while seed_segments_pcnt and (
            top_n_seeds < 1 or
            len(patterns_pcnt) < top_n_seeds
        ):
            (
                seed_pos,
                seed_segment
            ), pattern_pcnt = seed_segments_pcnt.most_common(1)[0]
            pattern_count = seed_segments_ct[(seed_pos, seed_segment)]
            # seed_segments.pop(seed_segment)

            selected_segments = [(seed_pos, seed_segment)]

            prev_segment = seed_segment
            for pos in range(
                seed_pos - self.segment_step,
                real_pos_start - 1,
                -self.segment_step
            ):
                pos_segments = segments_between_pcnt[pos]
                for segment, pcnt in pos_segments.most_common():
                    count = segments_between_ct[pos][segment]
                    if is_continuous(segment, prev_segment, self.segment_step):
                        selected_segments.append((pos, segment))
                        if pcnt < pattern_pcnt:
                            pattern_pcnt = pcnt
                        if count < pattern_count:
                            pattern_count = count
                        prev_segment = segment
                        # seed_segments.pop(segment, None)
                        # pos_segments.pop(segment)
                        break
                else:
                    break

            prev_segment = seed_segment
            for pos in range(
                seed_pos + self.segment_step,
                1 + real_pos_end,
                self.segment_step
            ):
                pos_segments = segments_between_pcnt[pos]
                for segment, pcnt in pos_segments.most_common():
                    count = segments_between_ct[pos][segment]
                    if is_continuous(prev_segment, segment, self.segment_step):
                        selected_segments.append((pos, segment))
                        if pcnt < pattern_pcnt:
                            pattern_pcnt = pcnt
                        if count < pattern_count:
                            pattern_count = count
                        prev_segment = segment
                        # seed_segments.pop(segment, None)
                        # pos_segments.pop(segment)
                        break
                else:
                    break

            node_map: Dict[Tuple[int, int], PosNA] = {}
            for segpos, segment in selected_segments:
                for node in segment:
                    if node is None:
                        continue
                    if node.pos < pos_start or node.pos > pos_end:
                        continue
                    nodekey = (node.pos, node.bp)
                    node_map[nodekey] = node
                seed_segments_pcnt[(segpos, segment)] -= pattern_pcnt
                segments_between_pcnt[segpos][segment] -= pattern_pcnt
                seed_segments_ct[(segpos, segment)] -= pattern_count
                segments_between_ct[segpos][segment] -= pattern_count
                if seed_segments_ct[(segpos, segment)] == 0:
                    del seed_segments_pcnt[(segpos, segment)]
                    del segments_between_pcnt[segpos][segment]
                    del seed_segments_ct[(segpos, segment)]
                    del segments_between_ct[segpos][segment]
            pattern_nodes = tuple(sorted(node_map.values()))
            if pattern_nodes and pattern_count is not None:
                patterns_pcnt[pattern_nodes] += pattern_pcnt
                patterns_ct[pattern_nodes] += pattern_count
        return {
            pattern: (patterns_ct[pattern], pcnt / 10000)
            for pattern, pcnt in patterns_pcnt.most_common()
        }

    @cython.ccall
    @cython.locals(count=cython.int)
    @cython.returns(cython.void)
    def add(
        self: 'SegFreq',
        segment: Tuple[Optional[PosNA], ...],
        count: int = 1
    ) -> None:
        """Add a segment to this object.

        :param segment: A tuple of nodes, where None represents a gap.
        :param count: The number of times to add the segment.
        """
        pos: int = cython.declare(cython.long, get_segment_pos(segment))
        if pos not in self._segments:
            self._segments[pos] = Counter()
        self._segments[pos][segment] += count
        self._max_segpos = self._max_segpos if self._max_segpos > pos else pos

    @cython.ccall
    @cython.returns(cython.void)
    def update(self: 'SegFreq', other: 'SegFreq') -> None:
        """Update this SegFreq with another SegFreq."""
        if self.segment_size != other.segment_size:
            raise ValueError("Incompatible segment sizes")
        for pos_segments in other._segments.values():
            for segment, count in pos_segments.items():
                self.add(segment, count)

    def dump(self: 'SegFreq', fp: TextIO) -> None:
        """Dump to a segfreq file."""
        csv_writer = csv.writer(fp)
        csv_writer.writerow([f'# segment_size={self.segment_size}'])
        csv_writer.writerow([f'# segment_step={self.segment_step}'])
        csv_writer.writerow([
            'pos', 'segment', 'offsets', 'count'])
        for pos, pos_segments in sorted(self._segments.items()):
            for segment, count in sorted(
                pos_segments.items(),
                key=lambda item: (
                    *(node.bp for node in item[0] if node),
                    *(node.na for node in item[0] if node)
                )
            ):
                csv_writer.writerow([
                    pos,
                    ''.join(chr(node.na) if node else '.' for node in segment),
                    ''.join(
                        '.' if
                        from_node is None or
                        to_node is None or
                        from_node.pos < to_node.pos
                        else '+'
                        for from_node, to_node in pairwise(segment)
                    ),
                    count
                ])

    @classmethod
    def load(cls, fp: TextIO) -> 'SegFreq':
        """Load from a segfreq file."""
        lines = fp.readlines()
        segment_size: int = DEFAULT_SEGMENT_SIZE
        segment_step: int = DEFAULT_SEGMENT_STEP
        nocmt_lines: Deque = deque()
        for line in lines:
            if line.startswith('#'):
                if re.match(r'^#\s*segment_size\s*=\s*\d+\s*$', line):
                    segment_size = int(line.split('=')[1].strip())
                elif re.match(r'^#\s*segment_step\s*=\s*\d+\s*$', line):
                    segment_step = int(line.split('=')[1].strip())
            else:
                nocmt_lines.append(line)

        csv_reader = csv.reader(nocmt_lines)
        next(csv_reader)
        self: 'SegFreq' = cls(segment_size, segment_step)
        for row in csv_reader:
            segment_list: List[Optional[PosNA]] = []
            prev_pos: int = int(row[0])
            prev_bp: int = 0
            for offset, na in zip('=' + row[2], row[1]):
                node: Optional[PosNA]
                if offset == '=':
                    if na == '.':
                        node = None
                    else:
                        node = PosNA(prev_pos, 0, ord(na))
                elif offset == '+':
                    node = PosNA(prev_pos, prev_bp + 1, ord(na))
                    prev_bp += 1
                else:
                    if na == '.':
                        node = None
                        prev_pos += 1
                        prev_bp = 0
                    else:
                        node = PosNA(prev_pos + 1, 0, ord(na))
                        prev_pos += 1
                        prev_bp = 0
                segment_list.append(node)
            self.add(tuple(segment_list), int(row[3]))
        return self
