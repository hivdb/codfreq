"""CLI group for working with CodFreq profile files."""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path
from typing import Annotated, Any, Iterable, cast

from Bio import SeqIO  # type: ignore[import-not-found]
from Bio.SeqFeature import CompoundLocation  # type: ignore[import-not-found]
import questionary  # type: ignore[import-not-found]
import typer
from pydantic import BaseModel, ConfigDict, ValidationError
from rich import print

from .codfreq_types import Profile


class GeneFeature(BaseModel):
    """Gene feature extracted from a GenBank record."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    ranges: list[tuple[int, int]]


class GenBankRecord(BaseModel):
    """Simplified representation of a GenBank record."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    accession: str
    sequence: str
    features: list[GeneFeature]


def _fetch_record(accession: str) -> GenBankRecord:  # pragma: no cover - network I/O
    """Retrieve a GenBank record for *accession* using Biopython.

    :param accession: GenBank accession identifier.
    :type accession: str
    :returns: Parsed GenBank record with sequence and gene features.
    :rtype: GenBankRecord
    :raises Exception: If the accession cannot be fetched.
    """

    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=nuccore&id={accession}&rettype=gb&retmode=text"
    )
    with urllib.request.urlopen(  # pragma: no cover - network I/O
        url
    ) as handle:  # type: ignore[no-untyped-call]
        record = SeqIO.read(handle, "genbank")  # type: ignore[no-untyped-call]

    features: list[GeneFeature] = []
    for feat in record.features:
        if feat.type not in {"gene", "CDS"}:
            continue
        gene_name = (
            feat.qualifiers.get("gene")
            or feat.qualifiers.get("product")
            or [None]
        )[0]
        if gene_name is None:
            continue
        location = feat.location
        parts: Iterable = (
            location.parts
            if isinstance(location, CompoundLocation)
            else [location]
        )
        ranges = [(int(p.start) + 1, int(p.end)) for p in parts]
        features.append(GeneFeature(name=cast(str, gene_name), ranges=ranges))

    return GenBankRecord(
        accession=accession, sequence=str(record.seq), features=features
    )


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
    default = [choice.value for choice in choices]
    selected: list[GeneFeature] = questionary.checkbox(
        "Select fragments to add:",
        choices=choices,
        default=default,  # type: ignore[arg-type]
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


def _auto_assembly_options(
    fragments: list[dict[str, object]]
) -> list[list[dict[str, object]]]:
    """Return possible assembly configurations derived from *fragments*.

    The function sorts fragments by genomic coordinates and fills gaps with
    inter-fragment regions. Overlaps are resolved by trimming either
    the left or right fragment. When overlaps exist, both trimming
    strategies are returned.

    ``trim`` entries are 1-based inclusive ranges relative to the fragment.

    :param fragments: Fragment configuration list.
    :type fragments: list[dict[str, object]]
    :returns: Candidate assembly configurations.
    :rtype: list[list[dict[str, object]]]
    """

    main = next((f for f in fragments if "refSequence" in f), None)
    if main is None:
        return []
    main_name = cast(str, main["fragmentName"])
    main_len = len(cast(str, main["refSequence"]))

    genes = [f for f in fragments if f.get("geneName")]
    if not genes:
        return [[{
            "name": main_name,
            "fromFragment": main_name,
            "refStart": 1,
            "refEnd": main_len,
        }]]

    def span(frag: dict[str, object]) -> tuple[int, int]:
        ranges = cast(list[tuple[int, int]], frag["refRanges"])
        starts = [r[0] for r in ranges]
        ends = [r[1] for r in ranges]
        return min(starts), max(ends)

    ordered = sorted(genes, key=lambda f: span(f)[0])

    def build(bias: str) -> list[dict[str, object]]:
        assemblies: list[dict[str, object]] = []
        prev_end = 0
        prev_gene = None
        prev_frag = None
        prev_idx = None

        for frag in ordered:
            start, end = span(frag)
            if start > prev_end + 1:
                assemblies.append(
                    {
                        "name": f"{prev_gene or 'start'}-{frag['geneName']}",
                        "fromFragment": main_name,
                        "refStart": prev_end + 1,
                        "refEnd": start - 1,
                    }
                )
            overlap = prev_end - start + 1 if prev_end >= start else 0
            entry: dict[str, object] = {"geneName": frag["geneName"]}
            if overlap > 0:
                if bias == "left":
                    entry["trim"] = [[1, overlap]]
                    start = prev_end + 1
                elif (
                    bias == "right"
                    and prev_idx is not None
                    and prev_frag is not None
                ):
                    prev_len = sum(
                        r[1] - r[0] + 1
                        for r in cast(
                            list[tuple[int, int]], prev_frag["refRanges"]
                        )
                    )
                    trim_range = [prev_len - overlap + 1, prev_len]
                    trim_list = cast(
                        list[Any], assemblies[prev_idx].setdefault("trim", [])
                    )
                    trim_list.append(trim_range)
                    prev_end = start - 1
            assemblies.append(entry)
            prev_end = end
            prev_gene = cast(str, frag["geneName"])
            prev_frag = frag
            prev_idx = len(assemblies) - 1

        if prev_end < main_len:
            assemblies.append(
                {
                    "name": f"{prev_gene}-end",
                    "fromFragment": main_name,
                    "refStart": prev_end + 1,
                    "refEnd": main_len,
                }
            )
        return assemblies

    left = build("left")
    right = build("right")
    options = [left]
    if right != left:
        options.append(right)
    return options


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


def validate_main() -> None:
    """CLI entry point for quick profile validation."""

    typer.run(validate)  # pragma: no cover - CLI entry point


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

    options = _auto_assembly_options(fragments)
    assemblies: list[dict[str, object]]
    if options:
        choice = options[0]
        if len(options) > 1:
            choice = questionary.select(
                "Select an assembly strategy:",
                choices=[
                    questionary.Choice(f"Option {i+1}", value=opt)
                    for i, opt in enumerate(options)
                ],
            ).ask()  # pragma: no cover - interactive selection
        print(f"Suggested assembly: {choice}")
        if questionary.confirm("Use this assembly configuration?").ask():
            assemblies = choice
        else:
            assemblies = _prompt_manual_assemblies()
    else:  # pragma: no cover - no fragments scenario
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
