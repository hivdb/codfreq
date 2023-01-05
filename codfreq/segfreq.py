import cython  # type: ignore
import re
import csv
from collections import deque

from typing import Dict, Tuple, Set, TextIO, List, Optional, Deque
from more_itertools import pairwise

from .posnas import PosNA

DEFAULT_SEGMENT_SIZE: int = 3


@cython.cfunc
@cython.inline
@cython.returns(cython.ulong)
def get_segment_pos(segment: Tuple[Optional[PosNA], ...]) -> int:
    for idx, node in enumerate(segment):
        if node is not None:
            return node.pos - idx
    return 0


@cython.cclass
class SegFreq:
    """A class to represent an NGS alignment in digraph form."""

    segment_size: int = cython.declare(cython.uint, visibility="public")
    pending_segments: Set[
        Tuple[Optional[PosNA], ...]
    ] = cython.declare(set, visibility="private")
    _nodes: Set[PosNA] = cython.declare(set, visibility="private")
    _edges: Dict[
        PosNA,
        Set[PosNA]
    ] = cython.declare(dict, visibility="private")
    _segments: Dict[
        int,
        Dict[Tuple[Optional[PosNA], ...], int]
    ] = cython.declare(dict, visibility="private")

    def __init__(self: 'SegFreq', segment_size: int):
        self.segment_size = segment_size
        self._nodes = set()
        self._edges = {}
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
    @cython.returns(set)
    def get_nodes(self: 'SegFreq') -> Set[PosNA]:
        """Return the nodes in the graph."""
        self.sync()
        return self._nodes

    nodes = property(get_nodes)

    @cython.ccall
    @cython.returns(dict)
    def get_edges(self: 'SegFreq') -> Dict[PosNA, Set[PosNA]]:
        """Return the edges in the graph."""
        self.sync()
        return self._edges

    edges = property(get_edges)

    @cython.ccall
    @cython.returns(cython.void)
    def sync(self: 'SegFreq') -> None:
        """Sync the graph with the pending segments."""
        if self.pending_segments:
            for segment in self.pending_segments:
                for node_from, node_to in pairwise(segment):
                    if node_from is None or node_to is None:
                        continue
                    self._nodes.add(node_from)
                    if node_from not in self._edges:
                        self._edges[node_from] = set()
                    self._edges[node_from].add(node_to)
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
            self._segments[pos] = {}
        if segment not in self._segments[pos]:
            self._segments[pos][segment] = 0
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
