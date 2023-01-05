import cython  # type: ignore
import csv

from typing import Dict, Tuple, Set, TextIO, Any, List
from more_itertools import pairwise


@cython.cclass
class PosNA:
    pos: int = cython.declare(cython.uint, visibility="public")
    bp: int = cython.declare(cython.uint, visibility="public")
    na: int = cython.declare(cython.uint, visibility="public")

    def __init__(self: 'PosNA', pos: int, bp: int, na: int):
        self.pos = pos
        self.bp = bp
        self.na = na

    def __hash__(self: 'PosNA') -> int:
        return hash((self.pos, self.bp, self.na))

    def __eq__(self: 'PosNA', other: Any) -> bool:
        if not isinstance(other, PosNA):
            return False
        return (self.pos, self.bp, self.na) == (other.pos, other.bp, other.na)

    def __repr__(self: 'PosNA') -> str:
        return f"PosNA({self.pos!r}, {self.bp!r}, {self.na!r})"


@cython.cclass
class SegFreq:

    nodes: Set[PosNA] = cython.declare(set, visibility="public")
    edges: Dict[
        PosNA,
        Set[PosNA]
    ] = cython.declare(dict, visibility="public")
    segments: Dict[
        Tuple[PosNA, ...],
        int
    ] = cython.declare(dict, visibility="public")

    def __init__(self: 'SegFreq'):
        self.nodes = set()
        self.edges = {}
        self.segments = {}

    def add(
        self: 'SegFreq',
        segment: Tuple[PosNA, ...],
        count: int = 1
    ) -> None:
        for node_from, node_to in pairwise(segment):
            self.nodes.add(node_from)
            if node_from not in self.edges:
                self.edges[node_from] = set()
            self.edges[node_from].add(node_to)

        if segment not in self.segments:
            self.segments[segment] = 0
        self.segments[segment] += count

    def update(self: 'SegFreq', other: 'SegFreq') -> None:
        for segment, count in other.segments.items():
            self.add(segment, count)

    def dump(self: 'SegFreq', fp: TextIO) -> None:
        csv_writer = csv.writer(fp)
        csv_writer.writerow([
            'pos', 'bp', 'segment', 'offsets', 'count'])
        for segment, count in sorted(
            self.segments.items(),
            key=lambda item: (
                *(node.pos for node in item[0]),
                *(node.bp for node in item[0]),
                *(node.na for node in item[0])
            )
        ):
            csv_writer.writerow([
                segment[0].pos,
                segment[0].bp,
                ''.join(chr(node.na) for node in segment),
                ''.join(
                    '+' if from_node.pos == to_node.pos else '.'
                    for from_node, to_node in pairwise(segment)
                ),
                count
            ])

    @classmethod
    def load(cls, fp: TextIO) -> 'SegFreq':
        csv_reader = csv.reader(fp)
        next(csv_reader)
        self: 'SegFreq' = cls()
        for row in csv_reader:
            segment_list: List[PosNA] = []
            start_node = PosNA(int(row[0]), int(row[1]), ord(row[2][0]))
            prev_node = start_node
            for offset, na in zip(row[3], row[2][1:]):
                if offset == '+':
                    node = PosNA(prev_node.pos, prev_node.bp + 1, ord(na))
                else:
                    node = PosNA(prev_node.pos + 1, 0, ord(na))
                segment_list.append(node)
                prev_node = node
            self.add(tuple(segment_list), int(row[4]))
        return self
