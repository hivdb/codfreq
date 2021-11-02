from typing import List, Optional, ByteString
from subprocess import Popen, PIPE


def compress(
    data: ByteString,
    compresslevel: int = 9,
    *,
    mtime: Optional[float] = None
) -> bytes:
    out: bytes
    error: bytes
    cmd: List[str] = ['pigz', f'-{compresslevel}', '-c']
    if mtime:
        cmd.extend(['-M', f'{mtime}'])
    proc: Popen = Popen(
        cmd,
        stdin=PIPE,
        stdout=PIPE,
        stderr=PIPE)
    out, error = proc.communicate(input=data)
    if error:
        error_text: str = error.decode('UTF-8').strip()
        if error_text:
            raise RuntimeError(error_text)
    return out


def decompress(data: ByteString) -> bytes:
    cmd: List[str] = ['pigz', '-d', '-c']
    proc: Popen = Popen(
        cmd,
        stdin=PIPE,
        stdout=PIPE,
        stderr=PIPE)
    out, error = proc.communicate(input=data)
    if error:
        error_text: str = error.decode('UTF-8').strip()
        if error_text:
            raise RuntimeError(error_text)
    return out
