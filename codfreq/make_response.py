import os
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Annotated
from collections.abc import Iterator

import typer


app = typer.Typer(pretty_exceptions_enable=False)


def utcnow_text() -> str:
    """Return the current UTC timestamp as an ISO formatted string."""
    return datetime.now(tz=timezone.utc).isoformat()


def yield_codfreqs(workdir: Path) -> Iterator[tuple[str, csv.DictReader]]:
    """Yield codfreq name and row iterator pairs from a directory.

    :param workdir: Directory containing CodFreq files.
    :returns: Generator of file name and CSV rows.
    """
    suffix = '.codfreq'
    for fname in os.listdir(workdir):
        if not fname.endswith(suffix):
            continue
        name = fname.rsplit(suffix, 1)[0]
        with open(os.path.join(workdir, fname), encoding='utf-8-sig') as fp:
            yield name, csv.DictReader(fp)


def yield_untrans(workdir: Path) -> Iterator[tuple[str, Any]]:
    """Yield untranslated region data from a directory.

    :param workdir: Directory containing untranslated region files.
    :returns: Generator of file name and JSON data.
    """
    suffix = '.untrans.json'
    for fname in os.listdir(workdir):
        if not fname.endswith(suffix):
            continue
        name = fname.rsplit(suffix, 1)[0]
        with open(os.path.join(workdir, fname), encoding='utf-8-sig') as fp:
            yield name, json.load(fp)


@app.command()
def make_response(
    workdir: Annotated[
        Path,
        typer.Argument(
            ..., exists=True, file_okay=False, dir_okay=True, resolve_path=True
        ),
    ],
    path_prefix: Annotated[str, typer.Argument(...)],
) -> None:
    """Create CodFreq file for response.

    :param workdir: Directory containing CodFreq outputs.
    :param path_prefix: Prefix used to derive the task key.
    :returns: None
    """
    uniqkey = path_prefix.split('/', 1)[-1]
    codfreqs: dict[str, list[dict[str, Any]]] = {}
    for name, rows in yield_codfreqs(workdir):
        name = f'{name}.codfreq'
        gpmap: dict[tuple[str, int], dict[str, Any]] = {}
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
    codfreq_list = [{
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
            'codfreqs': codfreq_list
        }))


if __name__ == '__main__':  # pragma: no cover
    app()  # pragma: no cover
