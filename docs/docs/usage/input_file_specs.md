---
title: Input file
layout: default
parent: Usage
nav_order: 1
permalink: docs/usage/input-file
---

# Input File

FluMut can analyze multiple Influenza A virus sequences simultaneously.
It can handle partial and complete genome sequences of multiple samples, of any HA (H1–H16) and NA (N1–N9) subtype.
You must provide one or more files containing the nucleotide sequences in FASTA format (an example can be downloaded [here](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/fasta_input_example.fa)).
Sequences must adhere to the [IUPAC code](https://www.bioinformatics.org/sms/iupac.html): gaps (`-`) are ignored, while any other character stops the analysis with an error.

## Samples and segments

FluMut relies on the FASTA header to assign each sequence to a sample.
For this reason, the header must contain a sample ID consistent among sequences of the same sample.
Sequences of the same sample can also be split across different FASTA files: they are merged as long as they share the same sample ID.

By default, the tool expects the sample ID followed by an underscore and the segment name (e.g. `my_sample_PB2`).

The segment does not need to be specified: FluMut aligns each sequence against the reference sequences in its database and picks the one that fits best.
The segment found in the header is used only as a hint to restrict the search, which makes the analysis faster:

- if it matches a segment name (`PB2`, `PB1`, `PA`, `HA`, `NP`, `NA`, `MP`, `NS`), only the references of that segment are considered
- if it matches a reference name (e.g. `H5`, `N1`), only that reference is considered
- if it matches neither, all references are considered

{: .note}
If the FASTA header does not match the regular expression, FluMut prints a warning and uses the whole header as sample ID.
Make sure your headers are parsed correctly, otherwise sequences of the same sample are reported as different samples.

## Custom FASTA header parsing

If your sequences have a different FASTA header composition, you can use the option `--name-regex` (or the `Name regex` field in the [GUI advanced options](usage-gui#advanced-options)).
This works with Python regular expressions, using either the named groups `sample` and `segment` or, if they are missing, a first capture group to find the sample and a second capture group to find the segment.
Here we provide some examples with most common cases:

| FASTA Header                 | Regular Expression                |
| ---------------------------- | --------------------------------- |
| my_sample_PB2                | `(.+)_(.+)`                       |
| my_sample\|PB2               | `(.+)\|(.+)`                      |
| my_sample-PB2                | `(.+)-(.+)`                       |
| my_sample-PB2_something_else | `(.+)-(.+?)_`                     |
| PB2_my_sample                | `(?P<segment>.+?)_(?P<sample>.+)` |
| my_sample (no segment)       | `(.+)`                            |

{: .note}
To find the regular expression that better fits your FASTA header you can try it on [Regex101](https://regex101.com/) selecting `Python` flavor.

## Quality checks

While translating the sequences, FluMut checks them for possible issues.
These are not errors: the analysis goes on, but the issues are reported as warnings and listed in the [`Checks` sheet](../output#checks) of the Excel output.

| Check | Description |
| -- | -- |
| Duplicate protein | the same protein is found more than once in a sample (e.g. two HA sequences with the same sample ID) |
| Frameshift | an insertion or deletion shifts the reading frame of a protein |
| Premature stop codon | a stop codon is found before the end of a protein |
| Elongation | the expected stop codon is missing at the end of a protein |
| No valid codon | the sequence does not cover any complete codon of a protein |
