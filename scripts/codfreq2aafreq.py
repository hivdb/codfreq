#! /usr/bin/env python3

import os
import csv
import gzip
import argparse
from collections import Counter, defaultdict

from typing import (
    Dict, TextIO, Tuple, Counter as tCounter,
    Optional, List, DefaultDict
)

CODON_TABLE: Dict[str, str] = {
    'TTT': 'F',
    'TTC': 'F',
    'TTA': 'L',
    'TTG': 'L',

    'CTT': 'L',
    'CTC': 'L',
    'CTA': 'L',
    'CTG': 'L',

    'ATT': 'I',
    'ATC': 'I',
    'ATA': 'I',
    'ATG': 'M',

    'GTT': 'V',
    'GTC': 'V',
    'GTA': 'V',
    'GTG': 'V',

    'TCT': 'S',
    'TCC': 'S',
    'TCA': 'S',
    'TCG': 'S',

    'CCT': 'P',
    'CCC': 'P',
    'CCA': 'P',
    'CCG': 'P',

    'ACT': 'T',
    'ACC': 'T',
    'ACA': 'T',
    'ACG': 'T',

    'GCT': 'A',
    'GCC': 'A',
    'GCA': 'A',
    'GCG': 'A',

    'TAT': 'Y',
    'TAC': 'Y',

    'CAT': 'H',
    'CAC': 'H',
    'CAA': 'Q',
    'CAG': 'Q',

    'AAT': 'N',
    'AAC': 'N',
    'AAA': 'K',
    'AAG': 'K',

    'GAT': 'D',
    'GAC': 'D',
    'GAA': 'E',
    'GAG': 'E',

    'TGT': 'C',
    'TGC': 'C',
    'TGG': 'W',

    'CGT': 'R',
    'CGC': 'R',
    'CGA': 'R',
    'CGG': 'R',

    'AGT': 'S',
    'AGC': 'S',
    'AGA': 'R',
    'AGG': 'R',

    'GGT': 'G',
    'GGC': 'G',
    'GGA': 'G',
    'GGG': 'G',

    'TAA': '*',
    'TGA': '*',
    'TAG': '*',
}


def codon_to_aa(codon: str) -> Tuple[str, str]:
    nogap = codon.replace('-', '')
    nogap_len = len(nogap)
    if codon in ('', '---'):
        return ('-', '')
    elif nogap_len >= 3:
        aa = CODON_TABLE[nogap[:3]]
        ins = ''
        for i in range(3, nogap_len, 3):
            triplet = nogap[i:i+3]
            if len(triplet) != 3:
                break
            ins += CODON_TABLE[triplet]
        return (aa, ins)
    else:
        return ('X', '')


def directory(path: str) -> str:
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f'{path} is not a directory')


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'codfreq_dir',
        type=directory,
        help='Directory containing CodFreq files'
    )
    parser.add_argument(
        'aafreq_dir',
        type=str,
        help='Directory to write AAFreq files'
    )
    args = parser.parse_args()
    os.makedirs(args.aafreq_dir, exist_ok=True)
    for filename in os.listdir(args.codfreq_dir):
        fpin: Optional[TextIO] = None
        fpout: Optional[TextIO] = None
        name: str
        lower = filename.lower()
        filepath = os.path.join(args.codfreq_dir, filename)
        try:
            if lower.endswith('.codfreq.gz'):
                name = filename[:-11]
                fpin = gzip.open(
                    filepath,
                    'rt',
                    encoding='utf-8-sig')
            elif lower.endswith('.codfreq'):
                name = filename[:-8]
                fpin = open(
                    filepath,
                    'r',
                    encoding='utf-8-sig')
            if fpin is None:
                continue
            aacounter: DefaultDict[
                Tuple[str, str],
                tCounter[Tuple[str, str]]
            ] = defaultdict(Counter)
            codons: DefaultDict[
                Tuple[str, str, str, str],
                List[str]
            ] = defaultdict(list)
            totals: Dict[Tuple[str, str], int] = {}
            for row in csv.reader(fpin):
                if not row:
                    continue
                elif row[0].startswith('#'):
                    continue
                elif row[0] == 'gene':
                    continue
                gene, pos, total, codon, count = row[:5]
                aa, ins = codon_to_aa(codon)
                aacounter[(gene, pos)][(aa, ins)] += int(count)
                codons[(gene, pos, aa, ins)].append(codon)
                totals[(gene, pos)] = int(total)

            fpout = open(
                os.path.join(args.aafreq_dir, name + '.aafreq.csv'),
                'w',
                encoding='utf-8-sig'
            )
            writer = csv.writer(fpout)
            writer.writerow([
                'gene', 'pos', 'total', 'aa', 'ins',
                'codons', 'count', 'pcnt'])
            for (gene, pos), counter in aacounter.items():
                totalint = totals[(gene, pos)]
                for (aa, ins), countint in counter.most_common():
                    codonstr = ';'.join(codons[(gene, pos, aa, ins)])
                    writer.writerow([
                        gene, pos, totalint, aa, ins,
                        codonstr, countint,
                        round(countint / totalint, 4)])
            print(fpout.name)
        finally:
            if fpin:
                fpin.close()
            if fpout:
                fpout.close()


if __name__ == '__main__':
    main()
