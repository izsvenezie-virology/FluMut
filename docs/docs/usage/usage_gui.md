---
title: GUI usage
layout: default
parent: Usage
nav_order: 2
permalink: docs/usage/usage-gui
---
TODO
# Basic Usage
You can get the output file in an [Excel format] (user-friendly) running:
```
flumut -x excel_output.xlsm your_fasta.fa
```
If you prefer the text outputs (machine-readable format) run:
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
