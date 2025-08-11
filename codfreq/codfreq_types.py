from typing import Literal, Any
from pydantic import BaseModel, ConfigDict, field_validator, model_validator

FASTQFileName = str
Header = str
SeqText = str
AAPos = int
NAPos = int
GeneText = str
NAChar = int
MultiNAText = bytes
CodonText = MultiNAText
AAChar = int
MultiAAText = bytes

NAPosRange = tuple[NAPos, NAPos]


class PairedFASTQ(BaseModel):
    """Metadata for a paired FASTQ sample.

    :param name: Sample identifier.
    :type name: Header
    :param pair: Tuple of R1/R2 FASTQ filenames. ``None`` for single-end reads.
    :type pair: tuple[FASTQFileName, FASTQFileName | None]
    :param n: Number of FASTQ files in the pair.
    :type n: int
    """

    model_config = ConfigDict(frozen=True, extra='forbid')

    name: Header
    pair: tuple[FASTQFileName, FASTQFileName | None]
    n: int


class Sequence(BaseModel):
    """FASTA sequence entry."""

    model_config = ConfigDict(frozen=True, extra='forbid')

    header: Header
    sequence: SeqText


class MainFragmentConfig(BaseModel):
    """Reference fragment configuration.

    :param fragmentName: Name of the fragment.
    :type fragmentName: Header
    :param refSequence: Reference nucleotide sequence.
    :type refSequence: SeqText
    """

    model_config = ConfigDict(frozen=True, extra='forbid')

    fragmentName: Header
    refSequence: SeqText


class CodonAlignmentConfig(BaseModel):
    """Configuration options for codon alignment.

    :param relRefStart: Start position relative to the reference.
    :type relRefStart: NAPos
    :param relRefEnd: End position relative to the reference.
    :type relRefEnd: NAPos
    :param windowSize: Sliding window size in amino acids.
    :type windowSize: AAPos | None
    :param minGapDistance: Minimum nucleotide distance between gaps.
    :type minGapDistance: NAPos | None
    :param relGapPlacementScore: Relative gap placement score string.
    :type relGapPlacementScore: str | None
    """

    model_config = ConfigDict(frozen=True, extra='forbid')

    relRefStart: NAPos
    relRefEnd: NAPos
    windowSize: AAPos | None = None
    minGapDistance: NAPos | None = None
    relGapPlacementScore: str | None = None


class DerivedFragmentConfig(BaseModel):
    """Fragment derived from a reference fragment.

    :param fragmentName: Name of the fragment.
    :type fragmentName: Header
    :param fromFragment: Source reference fragment name.
    :type fromFragment: Header
    :param geneName: Gene identifier if any.
    :type geneName: GeneText | None
    :param refRanges: Reference coordinate ranges.
    :type refRanges: list[NAPosRange]
    :param codonAlignment: Codon alignment configuration or ``False`` to
        disable alignment.
    :type codonAlignment: Literal[False] | list[CodonAlignmentConfig] | None
    """

    model_config = ConfigDict(frozen=True, extra='forbid')

    fragmentName: Header
    fromFragment: Header
    geneName: GeneText | None = None
    refRanges: list[NAPosRange]
    codonAlignment: Literal[False] | list[CodonAlignmentConfig] | None = None

    @model_validator(mode="before")
    @classmethod
    def merge_ref_ranges(cls, data: Any) -> Any:
        """Merge ``refStart``/``refEnd`` into ``refRanges`` if provided.

        :param data: Raw input data.
        :type data: Any
        :returns: Normalized data with ``refRanges`` populated.
        :rtype: Any
        """

        if isinstance(data, dict):
            refstart = data.pop("refStart", None)
            refend = data.pop("refEnd", None)
            if (
                "refRanges" not in data
                and isinstance(refstart, int)
                and isinstance(refend, int)
            ):
                data["refRanges"] = [(refstart, refend)]
        return data

    @field_validator("refRanges")
    @classmethod
    def ensure_ref_ranges(cls, value: list[NAPosRange]) -> list[NAPosRange]:
        """Ensure ``refRanges`` is not empty."""

        if not value:
            raise ValueError("refRanges cannot be empty")
        return value


FragmentConfig = MainFragmentConfig | DerivedFragmentConfig


class GeneAssemblyConfig(BaseModel):
    """Assembly configuration for a gene region.

    :param geneName: Gene identifier referenced in ``fragmentConfig``.
    :type geneName: str
    :param trim: Ranges to exclude from the assembled gene. Each tuple is a
        1-based inclusive interval relative to the fragment.
    :type trim: list[tuple[NAPos, NAPos]] | None
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    geneName: str
    trim: list[tuple[NAPos, NAPos]] | None = None

    @model_validator(mode="before")
    @classmethod
    def normalize_trim(cls, data: Any) -> Any:
        """Normalize ``trim`` entries to ``(start, end)`` tuples."""

        if (
            isinstance(data, dict)
            and "trim" in data
            and data["trim"] is not None
        ):
            norm = []
            for item in data["trim"]:
                if isinstance(item, int):
                    norm.append((item, item))
                else:
                    norm.append((item[0], item[1]))
            data["trim"] = norm
        return data


class RegionAssemblyConfig(BaseModel):
    """Inter-gene assembly configuration defined by explicit coordinates.

    :param name: Region name.
    :type name: str
    :param fromFragment: Source fragment name.
    :type fromFragment: str
    :param refStart: Start position in reference coordinates (1-based,
        inclusive).
    :type refStart: int
    :param refEnd: End position in reference coordinates (1-based,
        inclusive).
    :type refEnd: int
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    fromFragment: str
    refStart: int
    refEnd: int


SequenceAssemblyConfig = GeneAssemblyConfig | RegionAssemblyConfig


class NARegionConfig(BaseModel):
    """Definition of a nucleotide assembly region."""

    model_config = ConfigDict(frozen=True, extra='forbid')

    name: str
    fromFragment: str
    refStart: int
    refEnd: int


class RegionalConsensus(BaseModel):
    """Consensus sequence for an assembly region."""

    model_config = ConfigDict(frozen=True, extra='forbid')

    name: str
    refStart: NAPos
    refEnd: NAPos
    consensus: str


class Profile(BaseModel):
    """Top level profile configuration.

    :param version: Profile schema version.
    :type version: str
    :param fragmentConfig: Fragment configuration list.
    :type fragmentConfig: list[FragmentConfig]
    :param sequenceAssemblyConfig: Assembly region configuration list.
    :type sequenceAssemblyConfig: list[SequenceAssemblyConfig]
    """

    model_config = ConfigDict(frozen=True, extra='forbid')

    version: str
    fragmentConfig: list[FragmentConfig]
    sequenceAssemblyConfig: list[SequenceAssemblyConfig]

    @model_validator(mode="after")
    def check_assembly_continuity(self) -> "Profile":
        """Ensure ``sequenceAssemblyConfig`` covers the reference contiguously.

        The first region must start at position ``1`` of the main fragment and
        subsequent regions must abut without gaps or overlaps. The final region
        must end at the length of the main fragment.

        :raises ValueError: If the assembly does not cover the reference
            contiguously.
        """

        main = next(
            (
                f
                for f in self.fragmentConfig
                if isinstance(f, MainFragmentConfig)
            ),
            None,
        )
        if main is None or not self.sequenceAssemblyConfig:
            return self

        main_len = len(main.refSequence)

        gene_spans: dict[str, tuple[int, int]] = {}
        for frag in self.fragmentConfig:
            if isinstance(frag, DerivedFragmentConfig) and frag.geneName:
                starts = [r[0] for r in frag.refRanges]
                ends = [r[1] for r in frag.refRanges]
                gene_spans[frag.geneName] = (min(starts), max(ends))

        expected_start = 1
        for region in self.sequenceAssemblyConfig:
            if isinstance(region, GeneAssemblyConfig):
                if region.geneName not in gene_spans:
                    msg = (
                        f"Unknown gene '{region.geneName}' in "
                        "sequenceAssemblyConfig"
                    )  # pragma: no cover
                    raise ValueError(msg)  # pragma: no cover
                start, end = gene_spans[region.geneName]
            else:
                start, end = region.refStart, region.refEnd
            if start != expected_start:
                raise ValueError("sequenceAssemblyConfig is not continuous")
            expected_start = end + 1

        if expected_start - 1 != main_len:
            raise ValueError(
                "sequenceAssemblyConfig does not extend to end of reference"
            )  # pragma: no cover

        return self


class CodFreqRow(BaseModel):
    """Single CodFreq output row.

    :param gene: Gene name.
    :type gene: GeneText
    :param position: Amino acid position.
    :type position: AAPos
    :param total: Total codons observed at the position.
    :type total: int
    :param codon: Codon sequence encoded as bytes.
    :type codon: CodonText
    :param count: Number of reads supporting the codon.
    :type count: int
    :param total_quality_score: Sum of base qualities supporting the codon.
    :type total_quality_score: float
    """

    model_config = ConfigDict(frozen=True, extra='forbid')

    gene: GeneText
    position: AAPos
    total: int
    codon: CodonText
    count: int
    total_quality_score: float


#                                 refStart refEnd
#                                     v      v
FragmentInterval = tuple[list[tuple[NAPos, NAPos]], Header]

RefAAs = dict[AAPos, MultiAAText]
