import json
from typing import Any, Union


def dumps(obj: Any, *args: Any, **kwargs: Any) -> bytes:
    return json.dumps(obj, *args, **kwargs).encode()


def loads(b: Union[bytes, str], *args: Any, **kwargs: Any) -> Any:
    data = b.decode() if isinstance(b, bytes) else b
    return json.loads(data, *args, **kwargs)
