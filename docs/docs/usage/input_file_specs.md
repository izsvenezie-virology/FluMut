---
title: Input file
layout: default
parent: Usage
nav_order: 1
permalink: docs/usage/input-file
---

# Input File

FluMut can analyze multiple Influenza A virus sequences simultaneously.
It can handle partial and complete genome sequences of multiple samples.
You must provide a single file containing all the nucleotide sequences in FASTA format (an example can be downloaded [here](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/FASTA_input_example.fa)).
Sequences must adhere to the [IUPAC code](https://www.bioinformatics.org/sms/iupac.html).
FluMut relies on the FASTA header to assign the sequence to a specific sample.
For this reason, the header must contain a sample ID consistent among sequences of the same sample.
The segment can be autodetected or suggested in the FASTA header.

By default, the tool expects the sample ID followed by an underscore and the segment name (e.g. `my_sample_PB2`).

## Custom FASTA header parsing

If your sequences have a different FASTA header composition, you can use the option `--name-regex`.
This works with Python regular expressions, using a first capture group to find the sample and a second capture group to find the segment.
Here we provide some examples with most common cases:

| FASTA Header                 | Regular Expression                |
| ---------------------------- | --------------------------------- |
| my_sample_PB2                | `(.+)_(.+)`                       |
| my_sample\|PB2               | `(.+)\|(.+)`                      |
| my_sample-PB2                | `(.+)-(.+)`                       |
| my_sample-PB2_something_else | `(.+)-(.+?)_`                     |
| PB2_my_sample                | `(?P<segment>.+?)_(?P<sample>.+)` |

{: .note}
To find the regular expression that better fits your FASTA header you can try it on [Regex101](https://regex101.com/) selecting `Python` flavor.
