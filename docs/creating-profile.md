# Creating a Profile File

This document describes the structure of CodFreq profile files and how to author new ones.

## Profile schema

A profile is a JSON document validated against `codfreq.codfreq_types.Profile`. The top-level keys are:

- `version`: Schema version string.
- `fragmentConfig`: List of fragment definitions. Each element is either a **main fragment** or a **derived fragment**.
- `sequenceAssemblyConfig`: List of nucleotide or amino-acid assembly regions.

### MainFragmentConfig

Defines the base reference sequence.

```json
{
  "fragmentName": "Reference name",
  "refSequence": "ACGT..."
}
```

### DerivedFragmentConfig

Derived fragments reference a main fragment and specify coordinate ranges.

```json
{
  "fragmentName": "subfragment",
  "fromFragment": "Reference name",
  "refRanges": [[1, 100]],
  "geneName": "optional gene",
  "codonAlignment": [
    {
      "relRefStart": 1,
      "relRefEnd": 100,
      "minGapDistance": 15,
      "windowSize": null,
      "relGapPlacementScore": null
    }
  ]
}
```

### SequenceAssemblyConfig

Assembly regions describe how fragments are stitched together. Valid fields are:

- `name`
- `geneName`
- `fromFragment`
- `refStart`
- `refEnd`
- `trim`

The optional `trim` field records 1-based inclusive `(start, end)` ranges
relative to a fragment that should be excluded, typically to resolve overlaps.

All coordinate positions in the profile (for `refRanges`, `refStart`, `refEnd`,
and `trim`) are **1-based** and **inclusive**.

## Validating a profile

Run `validate-profile` to check a file:

```bash
validate-profile profiles/SARS2.json
```

The command prints validation errors and exits non-zero if the file does not match the schema.

To interactively build a profile, use:

```bash
profile create my_profile.json
```

You may provide a GenBank accession when prompted. The command downloads the
reference sequence with Biopython and displays all gene features in a checkbox
list that is preselected by default so you can deselect unwanted genes. Genes
with discontiguous ranges (e.g., the SEV glycoprotein) are represented as
multiple ``refRanges`` entries. After fragment selection, the tool suggests one
or more assembly configurations. Overlaps produce alternative strategies that
trim either the left or right gene; excluded ranges are captured in ``trim``.
Coordinate numbers such as ``refRanges`` pairs, ``refStart``, ``refEnd`` and
``trim`` are **1-based** and **inclusive**, matching the reference sequence
indexing used internally.

## Reference sequence checks

For fragments that embed a GenBank accession in the name (e.g. `Wuhan-Hu-1::NC_045512.2`), ensure the `refSequence` exactly matches the corresponding GenBank sequence and that all `refRanges` fall within the sequence length. Our cross-check found:

- `SARS2.json` and `SEV.json` match their GenBank references.
- `HIV1.json` uses a custom hybrid sequence and does not match the HXB2 accession.

## Steps to create a profile

1. Fetch the reference genome from GenBank.
2. Define a main fragment with the full reference sequence.
3. Add derived fragments with coordinate ranges for genes or regions of interest.
4. Add assembly regions to stitch fragments as needed.
5. Validate the JSON file with `validate-profile`.
