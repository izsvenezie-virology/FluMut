---
title: CLI usage
layout: default
parent: Usage
nav_order: 2
permalink: docs/usage/usage-cli
---

# CLI usage

## Basic Usage

You can get the output file in an [Excel format](../output#excel-output) (user-friendly) running:

```
flumut -x excel_output.xlsm your_fasta.fa
```

If you prefer the [text outputs](../output#markers-output) (machine-readable format) run:

```
flumut -m markers_output.tsv -M mutations_output.tsv -l literature_output.tsv your_fasta.fa
```

You can analyze several FASTA files at once, listing them all at the end of the command.
Sequences with the same sample ID are merged into the same sample, even if they come from different files (see [Input file](input-file#samples-and-segments)):

```
flumut -x excel_output.xlsm first.fa second.fa third.fa
```

At least one output must be requested, otherwise FluMut exits with an error without running the analysis.

## Update the database

FluMut uses a SQLite based database containing all information needed to perform the analysis.
The database contains the reference sequences, their annotation, the amino acid mutations and all the related marker information.

The database is shipped inside FluMut, and every new database is released as a new FluMut version.
You can check whether a newer FluMut is available with:

```
flumut --check-update
```

The recommended way to get the latest markers is to update FluMut itself (e.g. `pip install --upgrade flumut`).
If you cannot update FluMut, you can download the database of the latest release and use it with the `-D`/`--db-file` option:

```
flumut --download-db flumut_db.sqlite
flumut -D flumut_db.sqlite -x excel_output.xlsm your_fasta.fa
```

{: .warning}
A database can be used only with a FluMut version supporting its structure.
If the downloaded database is not compatible with your FluMut version, FluMut exits with an `IncompatibleVersionError`: in this case you need to update FluMut.

## Options

### General options

- `-h`/`--help`: print the help message and exit.
- `--version`: print FluMut version and exit.
- `--all-versions`: print FluMut and FluMutDB (database) versions and exit.
- `--check-update`: check on GitHub whether a newer FluMut is available and exit.
- `--download-db FILE`: download the database of the latest FluMut release to `FILE` and exit. See [Update the database](#update-the-database).
- `--loglevel [dbg|inf|wrn|err]`: set the verbosity of the messages printed on screen. Default is `wrn`, which reports warnings (such as the [quality checks](input-file#quality-checks)) and errors.

### Input options

- `-n`/`--name-regex`: change the regular expression used to parse FASTA headers.
  By default the regular expression used is `(?P<sample>.+)_(?P<segment>.+)`. More details can be found [here](input-file#custom-fasta-header-parsing).
- `-D`/`--db-file`: use a custom markers database instead of the bundled one.

### Analysis options

- `-r`/`--relaxed`: report all markers where at least one mutation is found. By default FluMut reports only markers if all mutations composing the markers are found in the input sequence.

### Output options

- `-m`/`--markers-output`: return as an output a text file containing the list of the identified markers (see [Markers output](../output#markers-output) for details).
- `-M`/`--mutations-output`: return as an output a text file containing the list of the identified mutations (see [Mutations output](../output#mutations-output) for details).
- `-l`/`--literature-output`: return as an output a text file containing the list of the references of the identified markers (see [Literature output](../output#literature-output) for details).
- `-x`/`--excel-output`: return as an output an excel file containing all the outputs, a summary table and the quality checks (see [Excel output](../output#excel-output) for details).
