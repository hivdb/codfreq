import cython  # type: ignore
import re
import csv
from collections import deque

from typing import Dict, Tuple, Set, TextIO, List, Optional, Deque
from more_itertools import pairwise

from .posnas import PosNA


@cython.cclass
class SegFreq:
    """A class to represent an NGS alignment in digraph form."""

    segment_size: int = cython.declare(cython.uint, visibility="public")
    segment_step: int = cython.declare(cython.uint, visibility="public")
    pending_segments: Set[
        Tuple[Optional[PosNA], ...]
    ] = cython.declare(set, visibility="private")
    _nodes: Set[PosNA] = cython.declare(set, visibility="private")
    _edges: Dict[
        PosNA,
        Set[PosNA]
    ] = cython.declare(dict, visibility="private")
    segments: Dict[
        Tuple[Optional[PosNA], ...],
        int
    ] = cython.declare(dict, visibility="public")

    def __init__(self: 'SegFreq', segment_size: int, segment_step: int):
        self.segment_size = segment_size
        self.segment_step = segment_step
        self._nodes = set()
        self._edges = {}
        self.segments = {}
        self.pending_segments = set()

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
        self.pending_segments.add(segment)

        if segment not in self.segments:
            self.segments[segment] = 0
        self.segments[segment] += count

    @cython.ccall
    @cython.returns(cython.void)
    def update(self: 'SegFreq', other: 'SegFreq') -> None:
        """Update the graph with another graph."""
        if self.segment_size != other.segment_size:
            raise ValueError("Incompatible segment sizes")
        if self.segment_step != other.segment_step:
            raise ValueError("Incompatible segment steps")
        for segment, count in other.segments.items():
            self.add(segment, count)

    def dump(self: 'SegFreq', fp: TextIO) -> None:
        """Dump the graph to a file."""
        fp.write(f'# segment_size={self.segment_size}\n')
        fp.write(f'# segment_step={self.segment_step}\n')
        csv_writer = csv.writer(fp)
        csv_writer.writerow([
            'pos', 'segment', 'offsets', 'count'])
        for segment, count in sorted(
            self.segments.items(),
            key=lambda item: (
                *(node.pos for node in item[0] if node),
                *(node.bp for node in item[0] if node),
                *(node.na for node in item[0] if node)
            )
        ):
            csv_writer.writerow([
                next(
                    node.pos - idx
                    for idx, node in enumerate(segment)
                    if node
                ),
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
        segment_size: int = 6
        segment_step: int = 1
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
