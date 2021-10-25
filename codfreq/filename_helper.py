import os
import re


def suggest_pair_name(fnpair, pattern):
    filename, _ = fnpair
    delimiter, offset, _, reverse = pattern
    dirpath, filename = os.path.split(filename)
    samfile = re.split(r'(?i)\.fastq(?:.gz)?$', filename)[0]
    if reverse == -1:
        return os.path.join(
            dirpath, samfile)
    samfile = samfile.split(delimiter)
    if reverse:
        samfile.reverse()
    samfile = samfile[:offset] + samfile[offset + 1:]
    if reverse:
        samfile.reverse()
    return os.path.join(
        dirpath, delimiter.join(samfile))


def name_file(fnpair, pattern, suffix):
    return suggest_pair_name(fnpair, pattern) + suffix


def name_bamfile(name, gene):
    return '{}.{}.bam'.format(name, gene)


def name_codfreq(name):
    return '{}.codfreq'.format(name)


def replace_ext(filename, toext, fromext=None, name_only=False):
    if name_only:
        filename = os.path.split(filename)[-1]
    if fromext:
        return filename[-len(fromext):] + toext
    else:
        return os.path.splitext(filename)[0] + toext
