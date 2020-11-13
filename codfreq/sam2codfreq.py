import os
import re
import csv
import click
from collections import defaultdict, Counter

from .paired_reads import iter_paired_reads
from .poscodons import iter_poscodons
from .codonutils import translate_codon
from .filename_helper import replace_ext
from .reference_helper import get_refaas

CODON_PATTERN = re.compile(r'^[ACGTYRWSKMBDHVN]{3}$')

CODFREQ_HEADER = [
    'gene', 'position',
    'total', 'codon', 'count',
    'refaa', 'aa', 'percent',
    'mean_quality_score'
]


def sam2codfreq(
    sampath,
    gene,
    refaas,
    reference_start=1,
    fnpair=None,
    log_format='text'
):
    codonfreqs = defaultdict(Counter)
    qualities = defaultdict(Counter)
    all_paired_reads = iter_paired_reads(sampath)
    for _, poscodons in iter_poscodons(all_paired_reads,
                                       reference_start=reference_start,
                                       fnpair=fnpair,
                                       description=sampath,
                                       log_format=log_format):
        for refpos, codon, qua in poscodons:
            codonfreqs[refpos][codon] += 1
            qualities[refpos][codon] += qua
        del poscodons
    for refpos, codons in sorted(codonfreqs.items()):
        total = sum(codons.values())
        for codon, count in codons.items():
            qua = qualities[refpos][codon] / count
            if CODON_PATTERN.match(codon):
                aa = translate_codon(codon)
            elif not codon or codon == '---':
                aa = 'del'
            elif len(codon) > 3:
                aa = 'ins'
            else:
                aa = 'X'
            yield {
                'gene': gene,
                'position': refpos,
                'total': total,
                'codon': codon,
                'refaa': refaas[refpos],
                'aa': aa,
                'count': count,
                'percent': count / total,
                'mean_quality_score': qua
            }


@click.command()
@click.argument(
    'workdir',
    type=click.Path(exists=True, file_okay=False,
                    dir_okay=True, resolve_path=True))
@click.option(
    '-g', '--gene',
    required=True,
    type=str,
    help='Gene name')
@click.option(
    '-r', '--reference',
    required=True,
    type=click.Path(exists=True, file_okay=True,
                    dir_okay=False, resolve_path=True))
@click.option(
    '--reference-start',
    type=int,
    default=1,
    show_default=True,
    help='Reference NA start (will be substract from the result)')
def main(workdir, gene, reference, reference_start):
    refaas = get_refaas(reference, reference_start)
    for dirpath, _, filenames in os.walk(workdir, followlinks=True):
        for fn in filenames:
            if fn[-4:].lower() not in ('.sam', '.bam'):
                continue
            sampath = os.path.join(dirpath, fn)
            codfreqfile = replace_ext(sampath, '.codfreq')
            with open(codfreqfile, 'w', encoding='utf-8-sig') as fp:
                writer = csv.DictWriter(fp, CODFREQ_HEADER)
                writer.writeheader()
                writer.writerows(sam2codfreq(
                    sampath,
                    gene=gene,
                    refaas=refaas,
                    reference_start=reference_start))


if __name__ == '__main__':
    main()
