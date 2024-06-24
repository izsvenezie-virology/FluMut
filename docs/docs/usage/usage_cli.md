---
title: CLI usage
layout: default
parent: Usage
nav_order: 1
permalink: docs/usage/usage-cli
---

# Basic Usage
You can get the output file in an [Excel format](../output#excel-output) (user-friendly) running:
```
flumut -x excel_output.xlsm your_fasta.fa
```
If you prefer the [text outputs](../output#markers-output) (machine-readable format) run:
```
flumut -m markers_output.tsv -M mutations_output.tsv -l literature_output.tsv your_fasta.fa
```

>**_Note_**: you should always [update](#update-marker-database) the database before any analysis to be sure to use the latest version available.

# Update marker database
FluMut uses a SQLite based database containing all information needed to perform the analysis.
The database contains the reference sequences, their annotation, the amino acid mutations and all the related marker information.

We will regularly update this database with new published markers.
You should always use the latest version of our database and you can do it just by running this command:
```
flumut --update
```


# Options

- `-v`/`--version`: print FluMut version.
- `--update`: update the database to the latest version.
- `--skip-unmatch-names`: when a FASTA header does not match the pattern of the regular expression FluMut skips this sequence and continues the analysis. By default, FluMut exits with an error.
- `--skip-unknown-segments`: when the segment name is not present in the database (eg. `>yoursample_P3`) FluMut skips the sequence and continues the analysis. By default, FluMut exits with an error.
- `-r`/`--relaxed`: report all markers where at least one mutation is found. By default FluMut reports only markers if all mutations composing the markers are found in the input sequence.
- `-n`/`--name-regex`: change the regular expression used to parse FASTA headers. 
    By default the regular expression used is (?P\<sample\>.+)_(?P\<segment\>.+). More details can be found [here](input-file#custom-fasta-header-parsing).
- `-D`/`--db-file`: use a custom markers database instead of the default one.
- `-m`/`--markers-output`: return as an output a text file containing the list of the identified markers (see [Markers output](../output#markers-output) for details).
- `-M`/`--mutations-output`: return as an output a text file containing the list of the identified mutations (see [Mutations output](../output#mutations-output) for details).
- `-l`/`--literature-output`: return as an output a text file containing the list of all the references in the database (see [Literature output](../output#literature-output) for details).
- `-x`/`--excel-output`: return as an output an excel file containing all the outputs and a summary table (see [Excel output](../output#excel-output) for details).
