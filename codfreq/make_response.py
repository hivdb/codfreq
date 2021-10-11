import os
import csv
import json
import click
from datetime import datetime, timezone


def utcnow_text():
    return datetime.now(tz=timezone.utc).isoformat()


def yield_codfreqs(workdir):
    suffix = '.codfreq'
    for fname in os.listdir(workdir):
        if not fname.endswith(suffix):
            continue
        name = fname.rsplit(suffix, 1)[0]
        with open(os.path.join(workdir, fname), encoding='utf-8-sig') as fp:
            yield name, csv.DictReader(fp)


def yield_untrans(workdir):
    suffix = '.untrans.json'
    for fname in os.listdir(workdir):
        if not fname.endswith(suffix):
            continue
        name = fname.rsplit(suffix, 1)[0]
        with open(os.path.join(workdir, fname), encoding='utf-8-sig') as fp:
            yield name, json.load(fp)


@click.command()
@click.argument(
    'workdir',
    type=click.Path(exists=True, file_okay=False,
                    dir_okay=True, resolve_path=True))
@click.argument('path_prefix', type=str)
def make_response(workdir, path_prefix):
    uniqkey = path_prefix.split('/', 1)[-1]
    """Create CodFreq file for response"""
    codfreqs = {}
    for name, rows in yield_codfreqs(workdir):
        name = '{}.codfreq'.format(name)
        gpmap = {}
        for row in rows:
            gene = row['gene']
            pos = int(row['position'])
            total = int(row['total'])
            codon = row['codon']
            if len(codon) < 3:
                continue
            count = int(row['count'])
            total_quality_score = float(row['total_quality_score'])
            if (gene, pos) not in gpmap:
                gpmap[(gene, pos)] = {
                    'gene': gene,
                    'position': pos,
                    'totalReads': total,
                    'allCodonReads': []
                }
            gpmap[(gene, pos)]['allCodonReads'].append({
                'codon': codon,
                'reads': count,
                'totalQualityScore': total_quality_score
            })
        codfreqs.setdefault(name, []).extend(gpmap.values())
    untrans_lookup = dict(yield_untrans(workdir))
    codfreqs = [{
        'name': name,
        'untranslatedRegions': untrans_lookup.get(
            name.rsplit('.codfreq', 1)[0]
        ),
        'allReads': all_reads
    } for name, all_reads in codfreqs.items()]
    with open(os.path.join(workdir, 'response.json'), 'w') as fp:
        fp.write(json.dumps({
            'taskKey': uniqkey,
            'lastUpdatedAt': utcnow_text(),
            'status': 'success',
            'codfreqs': codfreqs
        }))


if __name__ == '__main__':
    make_response()
