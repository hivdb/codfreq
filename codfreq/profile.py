"""CLI group for working with CodFreq profile files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Iterable, cast

from Bio import Entrez, SeqIO  # type: ignore[import-not-found]
from Bio.SeqFeature import CompoundLocation  # type: ignore[import-not-found]
import questionary  # type: ignore[import-not-found]
import typer
from pydantic import BaseModel, ConfigDict, ValidationError
from rich import print

from .codfreq_types import (
    Profile,
    FragmentConfig,
    DerivedFragmentConfig,
    MainFragmentConfig,
    GeneAssemblyConfig,
    RegionAssemblyConfig,
    SequenceAssemblyConfig,
)


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


def _fetch_record(
    accession: str, email: str
) -> GenBankRecord:
    """Retrieve a GenBank record for *accession* using Biopython's Entrez.

    :param accession: GenBank accession identifier.
    :type accession: str
    :param email: Email address required by NCBI for Entrez requests.
    :type email: str
    :returns: Parsed GenBank record with sequence and gene features.
    :rtype: GenBankRecord
    :raises Exception: If the accession cannot be fetched.
    """

    Entrez.email = email  # type: ignore[assignment]
    with Entrez.efetch(
        db="nuccore", id=accession, rettype="gb", retmode="text"
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
) -> list[FragmentConfig]:
    """Ask the user which gene features to keep as fragments.

    :param record: GenBank record containing features.
    :type record: GenBankRecord
    :param main_name: Name of the main fragment.
    :type main_name: str
    :returns: Derived fragments chosen by the user.
    :rtype: list[FragmentConfig]
    """

    if not record.features:
        return []
    choices = [
        questionary.Choice(
            f"{feat.name} {feat.ranges}", value=feat, checked=True
        )
        for feat in record.features
    ]
    selected: list[GeneFeature] = questionary.checkbox(
        "Select fragments to add:", choices=choices
    ).ask()
    fragments: list[FragmentConfig] = []
    for feat in selected:
        fragments.append(
            DerivedFragmentConfig(
                fragmentName=feat.name,
                fromFragment=main_name,
                geneName=feat.name,
                refRanges=feat.ranges,
            )
        )
    return fragments


def _prompt_manual_fragments(
    main_name: str | None,
) -> tuple[str, list[FragmentConfig]]:
    """Prompt the user for manually defined fragments.

    The first prompt collects the main fragment if *main_name* is ``None``. All
    subsequent fragments are treated as derived fragments referencing the main
    fragment. Coordinate ranges are **1-based** and **inclusive**.

    :param main_name: Existing main fragment name if already defined.
    :type main_name: str | None
    :returns: Tuple of the main fragment name and the list of entered
        fragments.
    :rtype: tuple[str, list[FragmentConfig]]
    """

    frags: list[FragmentConfig] = []
    if main_name is None:
        main_name = questionary.text("Main fragment name:").ask()
        seq = questionary.text(
            f"Reference sequence for {main_name}:"
        ).ask()
        frags.append(
            MainFragmentConfig(fragmentName=main_name, refSequence=seq)
        )

    while True:
        name = questionary.text(
            "Fragment name (leave blank to finish):"
        ).ask()
        if not name:
            break
        gene = questionary.text(
            "Gene name (same as fragment name if leave empty):"
        ).ask()
        ranges_text = questionary.text(
            "Reference ranges for this fragment (e.g., 1-5,8-10):"
        ).ask()
        ranges: list[tuple[int, int]] = []
        for part in ranges_text.split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                start_s, end_s = part.split("-", 1)
                ranges.append((int(start_s), int(end_s)))
            else:
                pos = int(part)
                ranges.append((pos, pos))
        frags.append(
            DerivedFragmentConfig(
                fragmentName=name,
                fromFragment=main_name,
                geneName=gene or name,
                refRanges=ranges,
            )
        )
    return main_name, frags


def _auto_assembly(
    fragments: list[FragmentConfig],
) -> list[SequenceAssemblyConfig]:
    """Generate an assembly configuration from *fragments*.

    Fragments are ordered by their genomic coordinates. Gaps between genes are
    filled by inter-gene regions, and overlaps are resolved by trimming the
    leftmost bases of the subsequent fragment. ``trim`` ranges are **1-based**
    and **inclusive**.

    :param fragments: Fragment configuration list.
    :type fragments: list[FragmentConfig]
    :returns: Assembly configuration using left-trimming for overlaps.
    :rtype: list[SequenceAssemblyConfig]
    """

    main = next(
        (f for f in fragments if isinstance(f, MainFragmentConfig)), None
    )
    if main is None:
        return []
    main_name = main.fragmentName
    main_len = len(main.refSequence)

    genes = [
        f
        for f in fragments
        if isinstance(f, DerivedFragmentConfig) and f.geneName
    ]
    if not genes:
        return [
            RegionAssemblyConfig(
                name=main_name,
                fromFragment=main_name,
                refStart=1,
                refEnd=main_len,
            )
        ]

    def span(frag: DerivedFragmentConfig) -> tuple[int, int]:
        starts = [r[0] for r in frag.refRanges]
        ends = [r[1] for r in frag.refRanges]
        return min(starts), max(ends)

    ordered = sorted(genes, key=lambda f: span(f)[0])

    assemblies: list[SequenceAssemblyConfig] = []
    prev_end = 0
    prev_gene: str | None = None
    for frag in ordered:
        start, end = span(frag)
        if start > prev_end + 1:
            assemblies.append(
                RegionAssemblyConfig(
                    name=f"{prev_gene or 'start'}-{frag.geneName}",
                    fromFragment=main_name,
                    refStart=prev_end + 1,
                    refEnd=start - 1,
                )
            )
        overlap = prev_end - start + 1 if prev_end >= start else 0
        trim = None
        if overlap > 0:
            trim = [(1, overlap)]
            start = prev_end + 1
        assemblies.append(
            GeneAssemblyConfig(geneName=cast(str, frag.geneName), trim=trim)
        )
        prev_end = end
        prev_gene = cast(str, frag.geneName)

    if prev_end < main_len:
        assemblies.append(
            RegionAssemblyConfig(
                name=f"{prev_gene}-end",
                fromFragment=main_name,
                refStart=prev_end + 1,
                refEnd=main_len,
            )
        )

    return assemblies


def _prompt_manual_assemblies(
    fragments: list[FragmentConfig],
) -> list[SequenceAssemblyConfig]:
    """Prompt the user to enter assembly regions manually.

    The user is asked whether each region is a gene or an inter-gene fragment.
    Inputs are validated to ensure regions are contiguous and cover the
    reference without gaps or overlaps.

    :param fragments: Available fragments for span lookup.
    :type fragments: list[FragmentConfig]
    :returns: User-provided assembly configuration.
    :rtype: list[SequenceAssemblyConfig]
    """

    gene_spans: dict[str, tuple[int, int]] = {}
    for frag in fragments:
        if isinstance(frag, DerivedFragmentConfig) and frag.geneName:
            starts = [r[0] for r in frag.refRanges]
            ends = [r[1] for r in frag.refRanges]
            gene_spans[frag.geneName] = (min(starts), max(ends))

    assemblies: list[SequenceAssemblyConfig] = []
    prev_end = 0
    while questionary.confirm("Add an assembly region?").ask():
        if questionary.confirm("Is this region a gene?").ask():
            gene_name = questionary.text("Gene name:").ask()
            if gene_name not in gene_spans:
                print(
                    f"[red]Unknown gene {gene_name}[/red]"
                )  # pragma: no cover - user input validation
                continue  # pragma: no cover - user input validation
            start, end = gene_spans[gene_name]
            if start != prev_end + 1:
                print(
                    "[red]Gap or overlap detected; re-enter region[/red]"
                )  # pragma: no cover - user input validation
                continue  # pragma: no cover - user input validation
            trim_text = questionary.text(
                "Trim ranges (e.g., 1-5,10)? leave blank for none:"
            ).ask()
            trim: list[tuple[int, int]] | None = None
            if trim_text:
                trim = []
                for part in trim_text.split(","):
                    part = part.strip()
                    if not part:
                        continue
                    if "-" in part:
                        a, b = part.split("-", 1)
                        trim.append((int(a), int(b)))
                    else:
                        pos = int(part)
                        trim.append((pos, pos))
            assemblies.append(
                GeneAssemblyConfig(geneName=gene_name, trim=trim)
            )
            prev_end = end
        else:
            region_name = questionary.text("Region name:").ask()
            from_fragment = questionary.text("Source fragment:").ask()
            ref_start = int(
                questionary.text(
                    "Reference start position (1-based, inclusive):",
                ).ask()
            )
            ref_end = int(
                questionary.text(
                    "Reference end position (1-based, inclusive):",
                ).ask()
            )
            if ref_start != prev_end + 1:
                print(
                    "[red]Gap or overlap detected; re-enter region[/red]"
                )  # pragma: no cover - user input validation
                continue  # pragma: no cover - user input validation
            assemblies.append(
                RegionAssemblyConfig(
                    name=region_name,
                    fromFragment=from_fragment,
                    refStart=ref_start,
                    refEnd=ref_end,
                )
            )
            prev_end = ref_end
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


validate_app = typer.Typer(
    pretty_exceptions_enable=False, help="Validate profile files"
)
validate_app.command()(validate)


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

    fragments: list[FragmentConfig] = []
    accession = questionary.text(
        "GenBank accession (leave blank for manual input):"
    ).ask()
    main_name: str | None = None
    if accession:
        try:
            email = questionary.text(
                "Email address for NCBI queries:",
            ).ask()
            record = _fetch_record(accession, email)
        except Exception as err:  # pragma: no cover - network failure
            print(f"[red]Failed to fetch {accession}: {err}[/red]")
            raise typer.Exit(code=1)
        main_name = questionary.text(
            "Main fragment name:", default=record.accession
        ).ask()
        fragments.append(
            MainFragmentConfig(
                fragmentName=main_name, refSequence=record.sequence
            )
        )
        fragments.extend(_prompt_genbank_fragments(record, main_name))
    main_name, manual_frags = _prompt_manual_fragments(main_name)
    fragments.extend(manual_frags)

    assemblies = _auto_assembly(fragments)
    if assemblies:
        print(f"Suggested assembly: {assemblies}")
        if not questionary.confirm(
            "Use this assembly configuration?"
        ).ask():
            assemblies = _prompt_manual_assemblies(
                fragments
            )
    else:
        assemblies = _prompt_manual_assemblies(fragments)

    profile_data = {
        "version": version,
        "fragmentConfig": [f.model_dump() for f in fragments],
        "sequenceAssemblyConfig": [a.model_dump() for a in assemblies],
    }

    try:
        Profile.model_validate(profile_data)
    except ValidationError as err:
        print("[red]Profile creation failed[/red]")
        for e in err.errors():
            loc = ".".join(str(item) for item in e["loc"])
            print(f"{loc}: {e['msg']}")
        raise typer.Exit(code=1)

    output.write_text(json.dumps(profile_data, indent=2), encoding="utf-8")
    print(f"[green]Profile written to {output}[/green]")


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    app()  # pragma: no cover - CLI entry point
