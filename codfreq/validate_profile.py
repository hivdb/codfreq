"""CLI for quick profile validation."""

from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer
from pydantic import ValidationError
from rich import print

from .codfreq_types import Profile


def main(
    profile: Annotated[
        Path,
        typer.Argument(
            ..., exists=True, file_okay=True, dir_okay=False, resolve_path=True
        ),
    ]
) -> None:
    """Validate a profile JSON file against the schema.

    :param profile: Path to the profile file.
    :type profile: Path
    :raises typer.Exit: If validation fails.
    """

    try:
        Profile.model_validate_json(profile.read_text(encoding="utf-8"))
    except ValidationError as err:
        print("[red]Profile validation failed[/red]")
        for e in err.errors():
            loc = ".".join(str(item) for item in e["loc"])
            print(f"{loc}: {e['msg']}")
        raise typer.Exit(code=1)
    print("[green]Profile is valid[/green]")


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    typer.run(main)  # pragma: no cover - CLI entry point
