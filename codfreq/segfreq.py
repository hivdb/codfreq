import cython  # type: ignore
import re
import csv
from collections import deque, Counter

from typing import (
    Dict, Tuple, Set, TextIO, List, Optional, Deque, Counter as tCounter
)

from more_itertools import pairwise

from .posnas import PosNA

DEFAULT_SEGMENT_SIZE: int = 3


@cython.cfunc
@cython.inline
@cython.returns(cython.ulong)
def get_segment_pos(segment: Tuple[Optional[PosNA], ...]) -> int:
    try:
        return segment[0].pos  # type: ignore
    except AttributeError:
        raise ValueError('Malformed segment: first node is None')


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def get_first_pos_nodes(
    segment: Tuple[Optional[PosNA], ...]
) -> Tuple[PosNA, ...]:
    first_pos: int = cython.declare(cython.ulong, get_segment_pos(segment))
    nodes: List[PosNA] = []
    try:
        for node in segment:
            if node.pos == first_pos:  # type: ignore
                nodes.append(node)  # type: ignore
            else:  # node.pos > first_pos
                break
    except AttributeError:
        raise ValueError('Malformed segment: first_pos node cannot be None')
    return tuple(nodes)


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def remove_first_pos(
    segment: Tuple[Optional[PosNA], ...]
) -> Tuple[Optional[PosNA], ...]:
    first_pos: int = cython.declare(cython.ulong, get_segment_pos(segment))
    for idx, node in enumerate(segment):
        if node is None or node.pos > first_pos:
            return segment[idx:]
    raise ValueError("Segment is malformed")


@cython.cfunc
@cython.inline
@cython.returns(tuple)
def remove_last_pos(
    segment: Tuple[Optional[PosNA], ...]
) -> Tuple[Optional[PosNA], ...]:
    if segment[-1] is None:
        return segment[:-1]
    else:
        last_pos = segment[-1].pos
        for idx, node in enumerate(reversed(segment)):
            if node is None or node.pos < last_pos:
                return segment[:len(segment) - idx]
        raise ValueError("Segment is malformed")


@cython.cfunc
@cython.inline
@cython.returns(cython.bint)
def is_continuous(
    left_segment: Tuple[Optional[PosNA], ...],
    right_segment: Tuple[Optional[PosNA], ...]
) -> bool:
    result: bool = (
        remove_first_pos(left_segment) ==
        remove_last_pos(right_segment)
    )
    return result


@cython.cclass
class SegFreq:
    """A class to represent an NGS alignment in digraph form."""

    segment_size: int = cython.declare(cython.uint, visibility="public")
    pending_segments: Set[
        Tuple[Optional[PosNA], ...]
    ] = cython.declare(set, visibility="private")
    _left_segments: Dict[
        Tuple[Optional[PosNA], ...],
        Tuple[Optional[PosNA], ...]
    ] = cython.declare(dict, visibility="private")
    _right_segments: Dict[
        Tuple[Optional[PosNA], ...],
        Set[Tuple[Optional[PosNA], ...]]
    ] = cython.declare(dict, visibility="private")
    _segments: Dict[
        int,
        tCounter[Tuple[Optional[PosNA], ...]]
    ] = cython.declare(dict, visibility="private")

    def __init__(self: 'SegFreq', segment_size: int):
        self.segment_size = segment_size
        self._left_segments = {}
        self._right_segments = {}
        self._segments = {}
        self.pending_segments = set()

    @cython.ccall
    @cython.locals(
        pos1=cython.ulong,
        pos2=cython.ulong,
        pos3=cython.ulong
    )
    @cython.returns(dict)
    def get_pos_codons(
        self: 'SegFreq',
        pos1: int,
        pos2: int = 0,
        pos3: int = 0
    ) -> Dict[bytes, int]:
        """Get the codon counts for a given position."""
        pos: int = cython.declare(cython.ulong)
        accessed: bool = cython.declare(cython.bint)
        codon: bytearray = cython.declare(bytearray)
        codon_bytes: bytes = cython.declare(bytes)
        if pos2 == 0:
            pos2 = pos1 + 1
        if pos3 == 0:
            pos3 = pos1 + 2
        if abs(pos2 - pos1) >= self.segment_size or \
                abs(pos3 - pos2) >= self.segment_size or \
                abs(pos3 - pos1) >= self.segment_size:
            raise ValueError('Positions are too far apart')
        codon_counts: Dict[bytes, int] = {}
        for segment, count in self._segments.get(pos1, {}).items():
            codon = bytearray()
            for pos in (pos1, pos2, pos3):
                accessed = False
                for node in segment:
                    if node and node.pos == pos:
                        codon.append(node.na)
                        accessed = True
                if not accessed:
                    break
            else:
                codon_bytes = bytes(codon)
                if codon_bytes not in codon_counts:
                    codon_counts[codon_bytes] = 0
                codon_counts[codon_bytes] += count
        return codon_counts

    @cython.ccall
    @cython.locals(
        pos_start=cython.ulong,
        pos_end=cython.ulong
    )
    @cython.returns(list)
    def get_consensus(
        self: 'SegFreq',
        pos_start: int,
        pos_end: int
    ) -> List[Optional[PosNA]]:
        """Get consensus segment for a range of positions."""
        pos: int = cython.declare(cython.ulong)
        consensus: List[Optional[PosNA]] = cython.declare(list, [])

        for pos in range(pos_start, pos_end + 1):
            if pos not in self._segments:
                consensus.append(None)
            else:
                pos_segments = self._segments[pos]
                ((segment, _), ) = pos_segments.most_common(1)
                consensus.extend([node for node in segment
                                  if node and node.pos == pos])
        return consensus

    @cython.ccall
    @cython.locals(
        pos_start=cython.ulong,
        pos_end=cython.ulong,
        min_count=cython.uint
    )
    @cython.returns(dict)
    def get_patterns_between(
        self: 'SegFreq',
        pos_start: int,
        pos_end: int,
        top_n: int
    ) -> Dict[Tuple[PosNA, ...], int]:
        """Get segments between two positions."""
        # self.sync()
        real_pos_end = pos_end - self.segment_size
        real_pos_end = pos_start if pos_start > real_pos_end else real_pos_end
        segments_between: Dict[
            int,
            tCounter[Tuple[Optional[PosNA], ...]]
        ] = {}
        seed_segments: tCounter[Tuple[Optional[PosNA], ...]] = Counter()
        for pos in range(pos_start, 1 + real_pos_end):
            segments_between[pos] = Counter(
                self._segments.get(pos, {}))
            seed_segments.update(segments_between[pos])

        patterns: Dict[Tuple[PosNA, ...], int] = {}
        pattern_count: int
        seed_segment: Tuple[Optional[PosNA], ...]
        prev_segment: Tuple[Optional[PosNA], ...]
        while seed_segments and len(patterns) < top_n:
            seed_segment, pattern_count = seed_segments.most_common(1)[0]
            # seed_segments.pop(seed_segment)

            selected_segments = [seed_segment]
            seed_pos = get_segment_pos(seed_segment)

            prev_segment = seed_segment
            for pos in range(seed_pos - 1, pos_start - 1, -1):
                pos_segments = segments_between[pos]
                for segment, count in pos_segments.most_common():
                    if is_continuous(segment, prev_segment):
                        selected_segments.append(segment)
                        if count < pattern_count:
                            pattern_count = count
                        prev_segment = segment
                        # seed_segments.pop(segment, None)
                        # pos_segments.pop(segment)
                        break
                else:
                    break

            prev_segment = seed_segment
            for pos in range(seed_pos, 1 + real_pos_end):
                pos_segments = segments_between[pos]
                for segment, count in pos_segments.most_common():
                    if is_continuous(prev_segment, segment):
                        selected_segments.append(segment)
                        if count < pattern_count:
                            pattern_count = count
                        prev_segment = segment
                        # seed_segments.pop(segment, None)
                        # pos_segments.pop(segment)
                        break
                else:
                    break

            node_map: Dict[Tuple[int, int], PosNA] = {}
            for segment in selected_segments:
                for node in segment:
                    if node is None:
                        continue
                    nodekey = (node.pos, node.bp)
                    node_map[nodekey] = node
                segpos = get_segment_pos(segment)
                seed_segments[segment] -= pattern_count
                segments_between[segpos][segment] -= pattern_count
                if seed_segments[segment] == 0:
                    del seed_segments[segment]
                    del segments_between[segpos][segment]
            pattern_nodes = tuple(sorted(node_map.values()))
            if pattern_count is not None:
                patterns[pattern_nodes] = pattern_count
        return patterns

    @cython.ccall
    @cython.returns(cython.void)
    def sync(self: 'SegFreq') -> None:
        """Sync the graph with the pending segments."""
        if self.pending_segments:
            for segment in self.pending_segments:
                left_seg = remove_first_pos(segment)
                right_seg = remove_last_pos(segment)
                self._left_segments[segment] = left_seg
                if right_seg not in self._right_segments:
                    self._right_segments[right_seg] = set()
                self._right_segments[right_seg].add(segment)
            self.pending_segments = set()

    @cython.ccall
    @cython.locals(count=cython.uint)
    @cython.returns(cython.void)
    def add(
        self: 'SegFreq',
        segment: Tuple[Optional[PosNA], ...],
        count: int = 1
    ) -> None:
        """Add a segment to the graph.

        :param segment: A tuple of nodes, where None represents a gap.
        :param count: The number of times to add the segment.
        """
        pos: int = cython.declare(cython.ulong, get_segment_pos(segment))
        self.pending_segments.add(segment)
        if pos not in self._segments:
            self._segments[pos] = Counter()
        self._segments[pos][segment] += count

    @cython.ccall
    @cython.returns(cython.void)
    def update(self: 'SegFreq', other: 'SegFreq') -> None:
        """Update the graph with another graph."""
        if self.segment_size != other.segment_size:
            raise ValueError("Incompatible segment sizes")
        for pos_segments in other._segments.values():
            for segment, count in pos_segments.items():
                self.add(segment, count)

    def dump(self: 'SegFreq', fp: TextIO) -> None:
        """Dump the graph to a file."""
        fp.write(f'# segment_size={self.segment_size}\r\n')
        csv_writer = csv.writer(fp)
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
        """Load a graph from a file."""
        lines = fp.readlines()
        segment_size: int = DEFAULT_SEGMENT_SIZE
        nocmt_lines: Deque = deque()
        for line in lines:
            if line.startswith('#'):
                if re.match(r'^#\s*segment_size\s*=\s*\d+\s*$', line):
                    segment_size = int(line.split('=')[1].strip())
            else:
                nocmt_lines.append(line)

        csv_reader = csv.reader(nocmt_lines)
        next(csv_reader)
        self: 'SegFreq' = cls(segment_size)
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
