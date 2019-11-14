#! /usr/bin/env python

from __future__ import print_function
from __future__ import division
import os
import re
import sys
import csv

from itertools import zip_longest
from collections import defaultdict

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

MIN_READ_DEPTH = int(os.environ.get('MIN_READ_DEPTH', 100))
MIN_ALLELE_COUNT = int(os.environ.get('MIN_ALLELE_COUNT', 5))


def main():
    if len(sys.argv) < 4:
        print("Usage: {} <PCNT_CUTOFF> "
              "<OUTPUT> <INPUT_DIRECTORY>".format(sys.argv[0]),
              file=sys.stderr)
        exit(1)
    pcnt_cutoff = float(sys.argv[1]) / 100
    with open(sys.argv[2], 'w') as out:
        inputdir = sys.argv[3]
        for fname in os.listdir(inputdir):
            lfname = fname.lower()
            if not lfname.endswith('.nucfreq') and \
                    not lfname.endswith('.nucfish'):
                continue
            fname = os.path.join(inputdir, fname)
            cons = []
            with open(fname, 'r') as fp:
                all_nas = defaultdict(set)
                reader = csv.reader(fp, delimiter='\t')
                for napos, total, na, read, insdetail, *_ in reader:
                    try:
                        napos = int(napos)
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
                    if na == 'i':
                        for ins in insdetail.split(', '):
                            insnas, insread = ins.split(' (')
                            insread = int(insread[:-1])
                            insfreq = insread / total
                            if insread < MIN_ALLELE_COUNT:
                                continue
                            if insfreq < pcnt_cutoff:
                                continue
                            all_nas[napos].add(insnas)
                    if na == 'd':
                        all_nas[napos].add('-')
                    else:
                        all_nas[napos].add(na)
                prev_pos = 0
                for pos, nas in sorted(all_nas.items()):
                    posgap = pos - prev_pos
                    if posgap > 1:
                        cons.append('N' * (posgap))
                    prev_pos = pos
                    for single_nas in zip_longest(*nas):
                        single_nas = set(single_nas)
                        single_nas -= {None, '-'}
                        single_nas = ''.join(sorted(single_nas))
                        single_nas = re.sub('[^ACGT]', '', single_nas)
                        if len(single_nas) > 1:
                            single_nas = REVERSED_AMBIGUOUS_NAS[single_nas]
                        cons.append(single_nas)
            cons = ''.join(cons).strip('N')
            header = fname.rsplit('/', 1)[-1].rsplit('.', 1)[0]
            out.write('>{}\n{}\n'.format(header, cons))


if __name__ == '__main__':
    main()
