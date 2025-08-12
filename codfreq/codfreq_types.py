from typing import Literal, Any, Annotated
from pydantic import (BaseModel, Field, ConfigDict,
                      field_validator, model_validator)

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

DEFAULT_CODON_ALIGN_WINDOW_SIZE = 10
DEFAULT_CODON_ALIGN_MIN_GAP_DISTANCE = 30


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
    :type windowSize: AAPos
    :param minGapDistance: Minimum nucleotide distance between gaps.
    :type minGapDistance: NAPos
    :param relGapPlacementScore: Relative gap placement score string.
    :type relGapPlacementScore: str | None
    """

    model_config = ConfigDict(frozen=True, extra='forbid')

    relRefStart: NAPos | None = None
    relRefEnd: NAPos | None = None
    windowSize: AAPos = DEFAULT_CODON_ALIGN_WINDOW_SIZE
    minGapDistance: NAPos = DEFAULT_CODON_ALIGN_MIN_GAP_DISTANCE
    relGapPlacementScore: str = ''


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
    :type codonAlignment: list[CodonAlignmentConfig] | Literal[False] | None
    """

    model_config = ConfigDict(frozen=True, extra='forbid')

    fragmentName: Header
    fromFragment: Header
    geneName: GeneText | None = None
    refRanges: list[NAPosRange]
    codonAlignment: list[CodonAlignmentConfig] | Literal[False] | None = None

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
    :type trim: list[tuple[NAPos, NAPos]]
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    geneName: str
    trim: Annotated[list[tuple[NAPos, NAPos]], Field(default_factory=list)]

    def __str__(self) -> str:
        """String representation of the gene assembly configuration."""
        trim_text = ','.join(
            f"-{start}-{end}" for start, end in self.trim
        ) or "all"
        return f"{self.geneName}:{trim_text}"

    @field_validator("trim", mode="before")
    @classmethod
    def normalize_trim(
        cls, value: list[tuple[int, int]] | list[int] | None
    ) -> list[tuple[int, int]] | None:
        """Normalize ``trim`` entries to ``(start, end)`` tuples."""

        if value is None:
            return None
        norm: list[tuple[int, int]] = []
        for item in value:
            if isinstance(item, int):
                norm.append((item, item))
            else:
                norm.append((item[0], item[1]))
        return norm


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

    def __str__(self) -> str:
        """String representation of the region."""
        return f"{self.name}[{self.fromFragment}]:{self.refStart}-{self.refEnd}"


SequenceAssemblyConfig = GeneAssemblyConfig | RegionAssemblyConfig


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

    version: Literal['20221213']
    fragmentConfig: list[
        Annotated[FragmentConfig, Field(union_mode='left_to_right')]
    ]
    sequenceAssemblyConfig: list[SequenceAssemblyConfig]

    @field_validator("fragmentConfig", mode="after")
    @classmethod
    def ensure_fragment_name_unique(
        cls, value: list[FragmentConfig]
    ) -> list[FragmentConfig]:
        """Ensure each fragment has a unique name."""
        names = set()
        for fragment in value:
            if fragment.fragmentName in names:
                raise ValueError(
                    f"Duplicate fragment name: {fragment.fragmentName}"
                )
            names.add(fragment.fragmentName)
        return value


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
