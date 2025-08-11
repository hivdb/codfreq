"""CLI group for working with CodFreq profile files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, cast

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import questionary
import typer
from pydantic import ValidationError
from rich import print

from .codfreq_types import Profile


def _fetch_record(accession: str) -> SeqRecord:
    """Retrieve a GenBank record for *accession*.

    :param accession: GenBank accession identifier.
    :type accession: str
    :returns: Parsed GenBank record.
    :rtype: SeqRecord
    :raises Exception: If the accession cannot be fetched.
    """

    Entrez.email = "anon@example.com"  # type: ignore[assignment]
    with Entrez.efetch(  # pragma: no cover - network I/O
        db="nuccore", id=accession, rettype="gb", retmode="text"
    ) as handle:  # type: ignore[no-untyped-call]
        return cast(
            SeqRecord, SeqIO.read(handle, "genbank")  # type: ignore[no-untyped-call]
        )


app = typer.Typer(pretty_exceptions_enable=False, help="Manage profile files")


@app.command()
def validate(
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


@app.command()
def create(
    output: Annotated[
        Path,
        typer.Argument(..., file_okay=True, dir_okay=False, resolve_path=True),
    ]
) -> None:
    """Interactively build a profile and write it to *output*.

    The prompts describe the expected fields. Coordinate positions such as
    ``refStart`` and ``refEnd`` are **1-based** and **inclusive**. If a
    GenBank accession is provided, the reference sequence and gene ranges are
    fetched automatically. Gene features spanning multiple ranges are
    represented using ``refRanges``.

    :param output: File to write the resulting profile JSON.
    :type output: Path
    :raises typer.Exit: If the generated profile fails validation.
    """

    version = questionary.text("Profile version:").ask()

    fragments: list[dict[str, object]] = []
    accession = questionary.text(
        "GenBank accession (leave blank for manual input):"
    ).ask()
    if accession:
        try:
            record = _fetch_record(accession)
        except Exception as err:  # pragma: no cover - network failure
            print(f"[red]Failed to fetch {accession}: {err}[/red]")
            raise typer.Exit(code=1)
        main_name = questionary.text(
            "Main fragment name:", default=str(record.id)
        ).ask()
        fragments.append(
            {"fragmentName": main_name, "refSequence": str(record.seq)}
        )
        for feat in record.features:
            if feat.type not in {"gene", "CDS"}:  # pragma: no cover - skip unrelated
                continue
            gene_name = (
                feat.qualifiers.get("gene")
                or feat.qualifiers.get("product")
                or [None]
            )[0]
            if gene_name is None:  # pragma: no cover - feature lacks gene name
                continue
            parts = getattr(feat.location, "parts", [feat.location])
            ranges: list[tuple[int, int]] = []
            for part in parts:
                start = int(part.start) + 1
                end = int(part.end)
                ranges.append((start, end))
            if questionary.confirm(
                f"Add derived fragment for {gene_name} with ranges {ranges}?"
            ).ask():
                fragments.append(
                    {
                        "fragmentName": gene_name,
                        "fromFragment": main_name,
                        "geneName": gene_name,
                        "refRanges": ranges,
                    }
                )

    while True:
        name = questionary.text(
            "Fragment name (leave blank to finish):"
        ).ask()
        if not name:
            break
        seq = questionary.text(
            f"Reference sequence for {name}:"
        ).ask()
        fragments.append({"fragmentName": name, "refSequence": seq})

    assemblies: list[dict[str, object]] = []
    while questionary.confirm("Add an assembly region?").ask():
        region_name = questionary.text("Region name:").ask()
        from_fragment = questionary.text("Source fragment:").ask()
        ref_start = int(
            questionary.text(
                "Reference start position (1-based, inclusive):"
            ).ask()
        )
        ref_end = int(
            questionary.text(
                "Reference end position (1-based, inclusive):"
            ).ask()
        )
        assemblies.append(
            {
                "name": region_name,
                "fromFragment": from_fragment,
                "refStart": ref_start,
                "refEnd": ref_end,
            }
        )

    profile_data = {
        "version": version,
        "fragmentConfig": fragments,
        "sequenceAssemblyConfig": assemblies,
    }

    try:
        Profile.model_validate(profile_data)
    except ValidationError as err:  # pragma: no cover - defensive validation
        print("[red]Profile creation failed[/red]")
        for e in err.errors():
            loc = ".".join(str(item) for item in e["loc"])
            print(f"{loc}: {e['msg']}")
        raise typer.Exit(code=1)

    output.write_text(json.dumps(profile_data, indent=2), encoding="utf-8")
    print(f"[green]Profile written to {output}[/green]")


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    app()  # pragma: no cover - CLI entry point
