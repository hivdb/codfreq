import cython  # type: ignore
import csv

from typing import Dict, Tuple, Set, TextIO, Any


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
class TriFreq:

    nodes: Set[PosNA] = cython.declare(set, visibility="public")
    edges: Dict[
        PosNA,
        Set[PosNA]
    ] = cython.declare(dict, visibility="public")
    triplets: Dict[
        Tuple[PosNA, PosNA, PosNA],
        int
    ] = cython.declare(dict, visibility="public")

    def __init__(self: 'TriFreq'):
        self.nodes = set()
        self.edges = {}
        self.triplets = {}

    def add_triplet(
        self: 'TriFreq',
        node_a: PosNA,
        node_b: PosNA,
        node_c: PosNA,
        count: int = 1
    ) -> None:
        self.nodes.add(node_a)
        self.nodes.add(node_b)
        self.nodes.add(node_c)

        if node_a not in self.edges:
            self.edges[node_a] = set()
        self.edges[node_a].add(node_b)

        if node_b not in self.edges:
            self.edges[node_b] = set()
        self.edges[node_b].add(node_c)

        triplet = (node_a, node_b, node_c)
        if triplet not in self.triplets:
            self.triplets[triplet] = 0
        self.triplets[triplet] += count

    def update(self: 'TriFreq', other: 'TriFreq') -> None:
        for triplet, count in other.triplets.items():
            self.add_triplet(*triplet, count)

    def dump(self: 'TriFreq', fp: TextIO) -> None:
        csv_writer = csv.writer(fp)
        csv_writer.writerow([
            'pos', 'triplet', 'bp1', 'bp2offset', 'bp3offset', 'count'])
        for (node_a, node_b, node_c), count in sorted(
            self.triplets.items(),
            key=lambda x: (
                x[0][0].pos,
                x[0][1].pos,
                x[0][2].pos,
                x[0][0].bp,
                x[0][1].bp,
                x[0][2].bp,
                x[0][0].na,
                x[0][1].na,
                x[0][2].na
            )
        ):
            csv_writer.writerow([
                node_a.pos,
                chr(node_a.na) + chr(node_b.na) + chr(node_c.na),
                node_a.bp,
                0 if node_b.pos == node_a.pos else '+',
                0 if node_c.pos == node_b.pos else '+',
                count
            ])

    @classmethod
    def load(cls, fp: TextIO) -> 'TriFreq':
        csv_reader = csv.reader(fp)
        next(csv_reader)
        self: 'TriFreq' = cls()
        for row in csv_reader:
            node_a = PosNA(int(row[0]), int(row[2]), ord(row[1][0]))
            if row[3] == '+':
                node_b = PosNA(node_a.pos + 1, 0, ord(row[1][1]))
            else:
                node_b = PosNA(node_a.pos, node_a.bp + 1, ord(row[1][1]))
            if row[4] == '+':
                node_c = PosNA(node_b.pos + 1, 0, ord(row[1][2]))
            else:
                node_c = PosNA(node_b.pos, node_b.bp + 1, ord(row[1][2]))
            self.add_triplet(node_a, node_b, node_c, int(row[5]))
        return self
