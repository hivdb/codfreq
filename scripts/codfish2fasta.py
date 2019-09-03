#! /usr/bin/env python

from __future__ import print_function
from __future__ import division
import re
import sys
import csv

from itertools import zip_longest
from collections import defaultdict

GENE_OFFSET = {
    'PR': -1,
    'PROTEASE': -1,
    'RT': 99 - 1,
    'IN': 99 + 560 - 1,
    'INT': 99 + 560 - 1
}

REVERSED_AMBIGUOUS_NAS = {
   'AT': 'W',
   'CG': 'S',
   'AC': 'M',
   'GT': 'K',
   'AG': 'R',
   'CT': 'Y',
   'CGT': 'B',
   'AGT': 'D',
   'ACT': 'H',
   'ACG': 'V',
   'ACGT': 'N',
}

MIN_READ_DEPTH = 100
MIN_ALLELE_COUNT = 5


def main():
    if len(sys.argv) < 4:
        print("Usage: {} <PCNT_CUTOFF> "
              "<OUTPUT> <IN1> <IN2> ...".format(sys.argv[0]),
              file=sys.stderr)
        exit(1)
    pcnt_cutoff = float(sys.argv[1]) / 100
    with open(sys.argv[2], 'w') as out:
        for fname in sys.argv[3:]:
            cons = []
            with open(fname, 'r') as fp:
                all_codons = defaultdict(set)
                reader = csv.reader(fp, delimiter='\t')
                for gene, cdpos, total, codon, read, *_ in reader:
                    try:
                        cdpos = int(cdpos)
                    except ValueError:
                        continue
                    read = int(read)
                    total = int(total)
                    if total == 0:
                        continue
                    freq = read / total
                    if total < MIN_READ_DEPTH:
                        continue
                    if read < MIN_ALLELE_COUNT:
                        continue
                    if freq < pcnt_cutoff:
                        continue
                    all_codons[cdpos + GENE_OFFSET[gene.upper()]].add(codon)
                prev_pos = 0
                for pos, codons in sorted(all_codons.items()):
                    posgap = pos - prev_pos
                    if posgap > 1:
                        cons.append('NNN' * (posgap))
                    prev_pos = pos
                    for nas in zip_longest(*codons):
                        nas = set(nas)
                        nas -= {None, '-'}
                        nas = ''.join(sorted(nas))
                        nas = re.sub('[^ACGT]', '', nas)
                        if len(nas) > 1:
                            nas = REVERSED_AMBIGUOUS_NAS[nas]
                        cons.append(nas)
            cons = ''.join(cons).strip('N')
            header = fname.rsplit('/', 1)[-1].rsplit('.', 1)[0]
            out.write('>{}\n{}\n'.format(header, cons))


if __name__ == '__main__':
    main()
