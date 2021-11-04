import os
import json
from typing import (
    Dict,
    List,
    Tuple,
    TextIO,
    BinaryIO,
    Optional,
    ByteString
)

from .cmdwrappers import pigz
from .codfreq_types import RegionalConsensus

import click  # type: ignore

EXT_UNTRANS_JSON = '.untrans.json'
EXT_CODFREQ = '.codfreq'


def find_codfreq_untrans_pairs(
    workdir: str
) -> List[Tuple[str, Optional[str]]]:
    key: str
    untrans_json: Optional[str]
    codfreq: Optional[str]
    pair: Tuple[Optional[str], Optional[str]]
    pairs: Dict[str, Tuple[Optional[str], Optional[str]]] = {}
    for dirpath, _, filenames in os.walk(workdir, followlinks=True):
        for filename in filenames:
            untrans_json = codfreq = None
            if filename.startswith('.'):
                continue
            elif filename.endswith(EXT_CODFREQ):
                key = filename[:-len(EXT_CODFREQ)]
                codfreq = os.path.join(dirpath, filename)
            elif filename.endswith(EXT_UNTRANS_JSON):
                key = filename[:-len(EXT_UNTRANS_JSON)]
                untrans_json = os.path.join(dirpath, filename)
            else:
                continue
            key = os.path.join(dirpath, key)
            if key in pairs:
                pair = pairs[key]
                if codfreq is None:
                    codfreq = pair[0]
                elif untrans_json is None:
                    untrans_json = pair[1]
            pairs[key] = (codfreq, untrans_json)
    return [
        (codfreq, untrans_json)
        for codfreq, untrans_json in pairs.values()
        if codfreq is not None
    ]


@click.command()
@click.argument(
    'workdir',
    type=click.Path(exists=True, file_okay=False,
                    dir_okay=True, resolve_path=True))
@click.option(
    '-v', '--verbose',
    is_flag=True,
    help='Verbose output')
def compress_codfreq(workdir: str, verbose: bool) -> None:
    codfreq: str
    untrans: Optional[str]
    fp: BinaryIO
    text_fp: TextIO
    untrans_objs: RegionalConsensus
    payload: ByteString
    pairs: List[Tuple[
        str, Optional[str]
    ]] = find_codfreq_untrans_pairs(workdir)
    for codfreq, untrans in pairs:
        payload = bytearray(b'\xef\xbb\xbf')  # UTF-8 bom
        if untrans:
            with open(untrans, 'rb') as fp:
                payload.extend(b'# --- untranslated regions begin ---\n')
                for untrans_obj in json.load(fp):
                    payload.extend(
                        '# {name} {refStart}..{refEnd}: {consensus}\n'
                        .format(**untrans_obj)
                        .encode('UTF-8')
                    )
                payload.extend(b'# --- untranslated regions end ---\n')
        with open(codfreq, 'r', encoding='UTF-8-sig') as text_fp:
            payload.extend(text_fp.read().encode('UTF-8'))
        payload = pigz.compress(payload)
        with open(codfreq + '.gz', 'wb') as fp:
            fp.write(payload)
        if verbose:
            click.echo(json.dumps({
                'op': 'compress-codfreq',
                'to': f'{codfreq}.gz'
            }))


if __name__ == '__main__':
    compress_codfreq()
