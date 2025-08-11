"""CLI group for working with CodFreq profile files."""

from __future__ import annotations

import json
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, cast

import questionary  # type: ignore[import-not-found]
import typer
from pydantic import ValidationError
from rich import print

from .codfreq_types import Profile


@dataclass
class GeneFeature:
    """A gene feature extracted from a GenBank record."""

    name: str
    ranges: list[tuple[int, int]]


@dataclass
class GenBankRecord:
    """Simplified representation of a GenBank record."""

    accession: str
    sequence: str
    features: list[GeneFeature]


def _download_genbank(accession: str) -> str:
    """Return the GenBank flat file for *accession*.

    :param accession: GenBank accession identifier.
    :type accession: str
    :returns: GenBank record in flat-file format.
    :rtype: str
    :raises Exception: If the accession cannot be fetched.
    """

    url = (  # pragma: no cover - network I/O
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype=gb&retmode=text"
    )
    with urllib.request.urlopen(  # pragma: no cover - network I/O
        url
    ) as handle:  # type: ignore[no-untyped-call]
        return cast(str, handle.read().decode("utf-8"))


def _parse_location(loc: str) -> list[tuple[int, int]]:
    """Parse a GenBank location string into numeric ranges."""

    loc = re.sub(r"complement\((.*)\)", r"\1", loc)
    if loc.startswith("join("):
        loc = loc[5:-1]
    ranges: list[tuple[int, int]] = []
    for part in loc.split(","):
        match = re.match(r"(\d+)\.\.(\d+)", part.strip())
        if match:
            ranges.append((int(match.group(1)), int(match.group(2))))
    return ranges


def _parse_genbank(text: str, accession: str) -> GenBankRecord:
    """Parse GenBank *text* into a :class:`GenBankRecord`.

    :param text: GenBank flat file content.
    :type text: str
    :param accession: Accession identifier.
    :type accession: str
    :returns: Parsed GenBank record with sequence and gene features.
    :rtype: GenBankRecord
    """

    feat_match = re.search(
        r"FEATURES\s+Location/Qualifiers\n(.*)\nORIGIN", text, re.DOTALL
    )
    feature_lines = feat_match.group(1).splitlines() if feat_match else []
    features: list[GeneFeature] = []
    i = 0
    while i < len(feature_lines):
        line = feature_lines[i]
        if line.startswith("     gene") or line.startswith("     CDS"):
            location = line[21:].strip()
            i += 1
            gene_name = None
            while i < len(feature_lines) and feature_lines[i].startswith(
                "                     "
            ):
                qual = feature_lines[i].strip()
                if qual.startswith("/gene="):
                    gene_name = qual.split("=", 1)[1].strip("\"")
                elif qual.startswith(
                    "/product="
                ) and gene_name is None:  # pragma: no cover
                    gene_name = qual.split("=", 1)[1].strip("\"")
                i += 1
            if gene_name:
                features.append(
                    GeneFeature(gene_name, _parse_location(location))
                )
            continue
        i += 1  # pragma: no cover - loop exit

    seq_match = re.search(r"ORIGIN\n(.*)\n//", text, re.DOTALL)
    seq_lines = seq_match.group(1).splitlines() if seq_match else []
    sequence = "".join(
        re.sub(r"[^acgtACGT]", "", ln) for ln in seq_lines
    ).upper()
    return GenBankRecord(accession, sequence, features)


def _fetch_record(accession: str) -> GenBankRecord:
    """Retrieve a GenBank record for *accession*.

    :param accession: GenBank accession identifier.
    :type accession: str
    :returns: Parsed GenBank record.
    :rtype: GenBankRecord
    :raises Exception: If the accession cannot be fetched.
    """

    text = _download_genbank(accession)
    return _parse_genbank(text, accession)


def _prompt_genbank_fragments(
    record: GenBankRecord, main_name: str
) -> list[dict[str, object]]:
    """Ask the user which gene features to keep as fragments.

    :param record: GenBank record containing features.
    :type record: GenBankRecord
    :param main_name: Name of the main fragment.
    :type main_name: str
    :returns: Derived fragments chosen by the user.
    :rtype: list[dict[str, object]]
    """

    if not record.features:  # pragma: no cover - no features
        return []
    choices = [
        questionary.Choice(f"{feat.name} {feat.ranges}", value=feat)
        for feat in record.features
    ]
    selected: list[GeneFeature] = questionary.checkbox(
        "Select fragments to add:", choices=choices
    ).ask()
    fragments: list[dict[str, object]] = []
    for feat in selected:
        fragments.append(
            {
                "fragmentName": feat.name,
                "fromFragment": main_name,
                "geneName": feat.name,
                "refRanges": feat.ranges,
            }
        )
    return fragments


def _prompt_manual_fragments() -> list[dict[str, object]]:
    """Prompt the user to enter additional fragments manually."""

    frags: list[dict[str, object]] = []
    while True:
        name = questionary.text(
            "Fragment name (leave blank to finish):"
        ).ask()
        if not name:
            break
        seq = questionary.text(f"Reference sequence for {name}:").ask()
        frags.append({"fragmentName": name, "refSequence": seq})
    return frags


def _auto_assembly(
    fragments: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Generate assembly regions from fragments with ``geneName``."""

    assemblies: list[dict[str, object]] = []
    for frag in fragments:
        gene = frag.get("geneName")
        if gene:
            assemblies.append({"geneName": gene})
    return assemblies


def _prompt_manual_assemblies() -> list[dict[str, object]]:
    """Prompt for assembly regions when auto generation is unsuitable."""

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
    return assemblies


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
            "Main fragment name:", default=record.accession
        ).ask()
        fragments.append(
            {"fragmentName": main_name, "refSequence": record.sequence}
        )
        fragments.extend(_prompt_genbank_fragments(record, main_name))

    fragments.extend(_prompt_manual_fragments())

    assemblies = _auto_assembly(fragments)
    print(f"Suggested assembly: {assemblies}")
    if not questionary.confirm("Use this assembly configuration?").ask():
        assemblies = _prompt_manual_assemblies()

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
