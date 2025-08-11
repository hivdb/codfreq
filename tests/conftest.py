from __future__ import annotations

from typing import Iterator

import pytest

from .mock_postalign import mock_postalign

_POSTALIGN = mock_postalign()
_POSTALIGN.__enter__()


@pytest.fixture(scope="session", autouse=True)
def _cleanup_postalign() -> Iterator[None]:
    """Undo the ``postalign``/``pysam``/``cython`` stubs after tests.

    Yields:
        Iterator[None]: ``None``.
    """
    yield
    _POSTALIGN.__exit__(None, None, None)
