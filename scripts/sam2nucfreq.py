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

# see https://en.wikipedia.org/wiki/Phred_quality_score
OVERALL_QUALITY_CUTOFF = int(os.environ.get('OVERALL_QUALITY_CUTOFF', 35))
LENGTH_CUTOFF = int(os.environ.get('LENGTH_CUTOFF', 50))
SITE_QUALITY_CUTOFF = int(os.environ.get('SITE_QUALITY_CUTOFF', 35))

NUM_PROCESSES = int(os.environ.get('NTHREADS', 2))
INPUT_QUEUE = Queue(NUM_PROCESSES)
OUTPUT_QUEUE = Queue()

PRODUCER_CAPACITY = 50000

CHUNKSIZE = 500

ERR_OK = 0b00
ERR_TOO_SHORT = 0b01
ERR_LOW_QUAL = 0b10


def get_na_counts(seq, qua, aligned_pairs, header):

    # pre-filter
    err = ERR_OK
    if len(seq) < LENGTH_CUTOFF:
        err |= ERR_TOO_SHORT
    if mean(qua) < OVERALL_QUALITY_CUTOFF:
        err |= ERR_LOW_QUAL
    if err:
        return None, err

    nas = defaultdict(list)
    prev_refpos = 0
    prev_seqpos0 = 0
    for seqpos0, refpos0 in aligned_pairs:

        if refpos0 is None:
            # insertion
            refpos = prev_refpos
        else:
            refpos = refpos0 + 1
            prev_refpos = refpos

        if seqpos0 is None:
            # deletion
            n, q = '-', qua[prev_seqpos0]
        else:
            n, q = seq[seqpos0], qua[seqpos0]
            prev_seqpos0 = seqpos0

        if refpos == 0:
            continue

        nas[refpos].append((n, q))

    result_nas = []
    for pos, nqs in sorted(nas.items()):
        qs = [q for _, q in nqs]
        meanq = sum(qs) / len(qs)
        if meanq < SITE_QUALITY_CUTOFF:
            continue
        result_nas.append((pos, ''.join(n for n, _ in nqs)))
    lastpos, lastna = result_nas[-1]
    # remove insertion at the end of sequence read
    result_nas[-1] = (lastpos, lastna[0])
    return result_nas, err


def reads_consumer():
    while True:
        chunk = INPUT_QUEUE.get()
        out_chunk = []
        for args in chunk:
            out_chunk.append(get_na_counts(*args))
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
                                       read.get_aligned_pairs(False))
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

    nafreqs = defaultdict(Counter)
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
                for refpos, na in results:
                    nafreqs[refpos][na] += 1

    with open(sys.argv[2], 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['Position', 'Total', 'NA', 'ReadCount'])
        for refpos, counter in sorted(nafreqs.items()):
            total = sum(counter.values())
            for na, read in sorted(counter.items(),
                                   key=lambda it: (-it[1], it[0])):
                writer.writerow([refpos, total, na, read])
    print('{} reads processed. Of them:'.format(num_finished))
    print('  Length of {} were too short'.format(num_tooshort))
    print('  Quality of {} were too low'.format(num_lowqual))
    print("{} done.".format(sys.argv[2]))


if __name__ == '__main__':
    main()
