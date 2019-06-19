#! /usr/bin/env python

from __future__ import print_function
import os
import sys
import csv
import pysam
from statistics import mean
from multiprocessing import Process, Queue
from collections import defaultdict, Counter

REF_CODON_OFFSET = 0

GENE_AA_RANGE = (
    ('PR', 56 + 1, 56 + 99),
    ('RT', 56 + 99 + 1, 56 + 99 + 560),
    ('IN', 56 + 99 + 560 + 1, 56 + 99 + 560 + 288)
)

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 25))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))
SITE_QUALITY_CUTOFF = int(os.environ.get('SITE_QUALITY_CUTOFF', 25))

NUM_PROCESSES = int(os.environ.get('NTHREADS', 2))
INPUT_QUEUE = Queue(NUM_PROCESSES)
OUTPUT_QUEUE = Queue()

PRODUCER_CAPACITY = 50000

CHUNKSIZE = 500

ERR_OK = 0b00
ERR_TOO_SHORT = 0b01
ERR_LOW_QUAL = 0b10


def get_codon_counts(seq, qua, aligned_pairs, header):

    # pre-filter
    err = ERR_OK
    if len(seq) < LENGTH_CUTOFF:
        err |= ERR_TOO_SHORT
    if mean(qua) < OVERALL_QUALITY_CUTOFF:
        err |= ERR_LOW_QUAL
    if err:
        return None, err

    codons = {}
    prev_seqpos = -1
    min_cdpos = aligned_pairs[-1][1] if aligned_pairs else 0xffffffff
    max_cdpos = 0
    for seqpos, refpos in aligned_pairs:
        if prev_seqpos > -1:
            codonpos = (refpos - 1 - REF_CODON_OFFSET) // 3 + 1
            codonbp = (refpos - 1 - REF_CODON_OFFSET) % 3
            min_cdpos = min(min_cdpos, codonpos)
            nas = []
            for n, q in zip(seq[prev_seqpos:seqpos],
                            qua[prev_seqpos:seqpos]):
                if q < SITE_QUALITY_CUTOFF:
                    n = '*'
                nas.append(n)
            if codonpos not in codons:
                codons[codonpos] = ['-'] * 3
            codons[codonpos][codonbp] = ''.join(nas)
        prev_seqpos = seqpos
    if prev_seqpos > -1:
        codonpos = (refpos - REF_CODON_OFFSET) // 3 + 1
        codonbp = (refpos - REF_CODON_OFFSET) % 3
        if codonpos not in codons:
            codons[codonpos] = ['-'] * 3
        codons[codonpos][codonbp] = seq[prev_seqpos]
        max_cdpos = codonpos
    results = []
    for cdpos in range(min_cdpos, max_cdpos + 1):
        codon = ''.join(codons.get(cdpos, ['---']))
        codon = codon if codon == '---' else codon.replace('-', '')
        # codon = codon.replace('-', '')
        if len(codon) - codon.count('*') < 3:
            continue
        codon = codon.replace('*', 'N')
        results.append((cdpos, codon))
    return results, err


def reads_consumer():
    while True:
        chunk = INPUT_QUEUE.get()
        out_chunk = []
        for args in chunk:
            out_chunk.append(get_codon_counts(*args))
        OUTPUT_QUEUE.put(out_chunk)


def reads_producer(filename, offset):
    with pysam.AlignmentFile(filename, 'rb') as samfile:
        chunk = []
        limit = offset + PRODUCER_CAPACITY
        for idx, read in enumerate(samfile.fetch()):
            if idx < offset:
                continue
            if idx >= limit:
                break
            seq, qua, aligned_pairs = (read.query_sequence,
                                       read.query_qualities,
                                       read.get_aligned_pairs(True))
            if len(chunk) == CHUNKSIZE:
                INPUT_QUEUE.put(chunk)
                chunk = []
            chunk.append((seq, qua, aligned_pairs, idx + offset))
        if chunk:
            INPUT_QUEUE.put(chunk)


def main():
    if len(sys.argv) != 3:
        print("Usage: {} <SAMFILE> <OUTPUT>".format(sys.argv[0]),
              file=sys.stderr)
        exit(1)
    with pysam.AlignmentFile(sys.argv[1], 'rb') as samfile:
        totalreads = samfile.count(until_eof=True)
    offset = 0
    while offset < totalreads:
        producer = Process(target=reads_producer, args=(sys.argv[1], offset))
        producer.daemon = True
        producer.start()
        offset += PRODUCER_CAPACITY

    for _ in range(NUM_PROCESSES):
        consumer = Process(target=reads_consumer)
        consumer.daemon = True
        consumer.start()

    codonfreqs = defaultdict(Counter)
    num_finished = 0
    num_tooshort = 0
    num_lowqual = 0
    while num_finished < totalreads:
        out_chunk = OUTPUT_QUEUE.get()
        for results, err in out_chunk:
            num_finished += 1
            if err & ERR_TOO_SHORT:
                num_tooshort += 1
            elif err & ERR_LOW_QUAL:
                num_lowqual += 1
            else:
                for cdpos, codon in results:
                    codonfreqs[cdpos][codon] += 1

    with open(sys.argv[2], 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        for cdpos, counter in sorted(codonfreqs.items()):
            for gene, start, end in GENE_AA_RANGE:
                if start <= cdpos <= end:
                    cdpos -= start - 1
                    break
            else:
                continue
            total = sum(counter.values())
            for codon, read in sorted(counter.items(),
                                      key=lambda it: (-it[1], it[0])):
                writer.writerow([gene, cdpos, total, codon, read])
    print('{} reads processed. Of them:'.format(num_finished))
    print('  Length of {} were too short'.format(num_tooshort))
    print('  Quality of {} were too low'.format(num_lowqual))
    print("{} done.".format(sys.argv[2]))


if __name__ == '__main__':
    main()
