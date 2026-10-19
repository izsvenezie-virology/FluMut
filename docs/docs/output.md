---
title: Output
layout: default
nav_order: 4
permalink: docs/output
---

# Output

## Excel output

This is the most user-friendly and complete output.
You can obtain this output using the [`--excel-output` option](usage/usage-cli#options).
You can download the example file [here](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/excel_output_example.xlsm).

It contains 6 sheets:

- [Checks](#checks)
- [Markers per sample](#markers-per-sample)
- [Samples per marker](#samples-per-marker)
- [Mutations](#mutations)
- [Markers](#markers)
- [Literature](#literature)

{: .note}
To enable the navigation feature the extension of the Excel file must be `.xlsm`.
If you don't care about navigation, you can use `.xlsx` extension.

### Checks

`Checks` reports the FluMut and database versions used for the analysis, followed by the list of issues found by the [quality checks](usage/input-file#quality-checks) (e.g. frameshifts, premature stop codons, duplicated proteins).
If the list is empty no issue was detected in your sequences.

{: .important}
Always check this sheet: an issue in a sequence may affect the reliability of the markers identified in that sample.

{: .note}
This sheet also reports the FluMut and database version used to generate it.


### Markers per sample

`Markers per sample` reports the list of the identified markers grouped by sample.
Specifically, it contains six columns:

- **_Sample_**: sample ID;
- **_Marker_**: markers identified in the sequences.
  If the `--relaxed` option is selected, markers for which at least a mutation composing the marker is identified are reported.
  Markers flagged with `^` are [high priority markers](#high-priority-markers).
- **_Mutations in your sample_**: mutations composing the marker that are identified in the sequences analysed.
  This column is identical to the column `Marker` unless the [`--relaxed` option](usage/usage-cli#options) was selected.
  With the `--relaxed` option this column allows to inspect which of the mutations composing the marker are identified in the analysed sequences.
- **_Effect_**: biological effects of the marker, as described in the literature.
- **_Subtype_**: influenza subtypes where the effect of the marker has been observed or tested.
- **_Literature_**: bibliographic references where the marker has been described.

{: .note}
By double-clicking on cells in the `Marker` or `Mutations in your sample` columns, you are redirected to the `Mutations` sheet, with columns filtered to display only mutations of selected marker.

{: .note}
By double-clicking on cells containing references in the `Literature` column you are redirected to the `Literature` sheet.
Rows are filtered to display only papers describing the selected effect.

![](../images/xls-mps.jpg)

### Samples per marker

This sheet displays the number of samples in which the marker is observed (column `Samples count`).
If [`--relaxed` option](usage/usage-cli#options) was selected, samples for which at least one mutation of the marker are counted.

{: .note}
By double-clicking on cells in the `Marker` or `Mutations in your sample` columns, you are redirected to the `Mutations` sheet, with columns filtered to display only mutations of selected marker.

{: .note}
By double-clicking on cells containing references in the `Literature` column you are redirected to the `Literature` sheet.
Rows are filtered to display only papers describing the selected effect.

![](../images/xls-spm.jpg)

### Mutations

This sheet shows the amino acid present in the positions of mutations of interest for each sample.
Only mutations identified in at least one of the analysed samples are reported.
The mutations are ordered by segment, then by protein, then by position and then by amino acid mutation.

| blank cell | no sequence available for the protein |
| `-` | deletion |
| `?` | undefined amino acid (at least one N nucleotide in the triplet, or incomplete codon) |
| `*` | stop codon |
| two or more amino acids | codon codes for alternative amino acids (nucleotide degeneration within the codon causing a nonsynonymous mutation) |

![](../images/xls-muts.jpg)

### Markers

The `Markers` sheet lists for each sample, the marker and mutation identified.
If the marker has multiple effects or it was tested on multiple subtypes it will be composed of different rows (one row per effect/subtype).
It contains the same columns described in [Markers per sample](#markers-per-sample).

![](../images/xls-mark.jpg)

### Literature

The `Literature` sheet reports the references describing the markers identified in the analysed samples, with title, authors, year, journal, link and DOI.

![](../images/xls-lit.jpg)

## High priority markers

Some markers are flagged as high priority: these are the markers listed in Table 4 of the [EFSA and ECDC (2025)](https://doi.org/10.2903/j.efsa.2025.9191) scientific report.
In the Excel and in the [Markers output](#markers-output) their name is followed by `^` (e.g. `PB2:627K^`), so that you can spot them at a glance.

## Markers output

This output is obtained with the [`--markers-output` option](usage/usage-cli#options).
The file reports in a tab-delimited value format all the information present in [Markers](#markers) sheet of the excel output, described above.

## Mutations output

This output is obtained with the [`--mutations-output` option](usage/usage-cli#options).
The file reports in a tab-delimited value format all the information present in [Mutations](#mutations) sheet of the excel output, described above.

## Literature output

This output is obtained with the [`--literature-output` option](usage/usage-cli#options).
The file reports in a tab-delimited value format all the information present in [Literature](#literature) sheet of the excel output, described above.

{: .note}
The quality checks are reported only in the [Checks](#checks) sheet of the Excel output.
When using text outputs, they are printed on screen as warnings.
