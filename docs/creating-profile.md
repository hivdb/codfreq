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

Existing profiles include an additional `trim` field which is **not** part of the current schema and causes validation to fail.

All coordinate positions in the profile (for `refRanges`, `refStart`, and
`refEnd`) are **1-based** and **inclusive**.

## Validating a profile

Run `profile validate` to check a file:

```bash
profile validate profiles/SARS2.json
```

The command prints validation errors and exits non-zero if the file does not match the schema.

To interactively build a profile, use:

```bash
profile create my_profile.json
```

You may provide a GenBank accession when prompted. The command downloads the
reference sequence and suggests derived fragments for each gene feature. Genes
with discontiguous ranges (e.g., the SEV glycoprotein) are represented as
multiple ``refRanges`` entries. Coordinate numbers such as ``refRanges`` pairs,
``refStart``, and ``refEnd`` are **1-based** and **inclusive**, matching the
reference sequence indexing used internally.

## Reference sequence checks

For fragments that embed a GenBank accession in the name (e.g. `Wuhan-Hu-1::NC_045512.2`), ensure the `refSequence` exactly matches the corresponding GenBank sequence and that all `refRanges` fall within the sequence length. Our cross-check found:

- `SARS2.json` and `SEV.json` match their GenBank references.
- `HIV1.json` uses a custom hybrid sequence and does not match the HXB2 accession.

## Steps to create a profile

1. Fetch the reference genome from GenBank.
2. Define a main fragment with the full reference sequence.
3. Add derived fragments with coordinate ranges for genes or regions of interest.
4. Add assembly regions to stitch fragments as needed.
5. Validate the JSON file with `profile validate`.
