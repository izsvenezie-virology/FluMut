# MutFinder
[![Latest Release](https://gitlab.izsvenezie.it/biv/flumut/mutfinder/-/badges/release.svg)](https://gitlab.izsvenezie.it/biv/flumut/mutfinder/-/releases)

A tool to easly find mutations of interest in influenza viruses.

## Getting started

### Prerequisites
Before installing MutFinder you must have:
- [Python](https://www.python.org/downloads/)
- [Pip](https://pypi.org/project/pip/)

### Installation
- Download latest release zip from [MutFinder repository](https://gitlab.izsvenezie.it/biv/flumut/mutfinder/-/releases)
- Install using pip:
    ```
    pip install mutfinder.zip
    ```
- Check installation:
    ```
    mutfinder --help
    ```

## Usage
### Basic usage
```
mutfinder -x excel_output.xlsm -t markers_output.tsv -m mutations_output.tsv your_fasta.fa
```
The above command parse ```your_fasta.fa``` file and produces all possible [outputs](#outputs).

### Update database
```
mutfinder --update
```
This command updates to latest version of marker's database from local NAS.
In the future this command will download the database from a repository.

### Try examples
We provide two small datasets to try the program.

```
unzip mutfinder.zip
cd mutfinder/examples
mutfinder -t markers_output.tsv -m mutations_output.tsv single_sample.fa
mutfinder -x excel_output.xlsm multiple_samples.fa 
```
These commands produces 2 text outputs for single_sample.fa and the Excel output for multiple_samples.fa.

## Input
MutFinder takes a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file as input.
Sequences must be nucleotidic and must contains only [IUPAC nucleotide codes](https://www.bioinformatics.org/sms/iupac.html).
The file should contains all segments for each sample in order to find markers with mutations on different segments.
MutFinder is capable to analyze multiple samples at once.

If you have your samples or your segments splitted on multiple files you can perform analysis with something like:
```
cat *.fa | mutfinder -x excel_output.xlsm -
```

## Outputs
MutFinder can produce 3 different outputs:
- [Excel](#excel)
- [Tabular](#tabular)
- [Matrix](#matrix)

### Excel
Options: ```-x FILENAME```/```--excel-output FILENAME```
```
mutfinder -x excel_output.xlsm your_fasta.fa
mutfinder --excel-output excel_output.xlsm your_fasta.fa
```

This is the most user-friendly output and contains both [Tabular](#tabular) and [Matrix](#matrix) outputs.
The output file contains 5 sheets:
- Markers per sample: report of markers for each sample
- Samples per marker: report sample count for each marker
- Mutations: contains [Matrix](#matrix) output
- Markers: contains [Tabular](#tabular) output
- Papers: contains the list of papers from the database

> **_NOTE:_** in "Markers per sample" and "Samples per marker" by double clicking a cell in "Markers mutations" or "Found mutations" columns you will be redirected to "Mutations" sheet.
> Columns are filtered to display only mutations in the selected marker.
>
> By double clicking a cell in "Papers" column you will be redirected to "Papers" sheet. 
> Rows are filtered to display only papers in selected effect.

### Tabular
Options: ```-t FILENAME```/```--tabular-output FILENAME```
```
mutfinder -t markers_output.tsv your_fasta.fa
mutfinder --tabular-output markers_output.tsv your_fasta.fa
```

Outputs a tab-separeted file with the following columns:
- Sample: the sample ID
- Marker mutations: list of mutations making up the marker
- Found mutations: list of mutations making up the marker actually found in the sample
- Effect: reported effect for the specific marker
- Subtype: subtype on which the effect was studied
- Papers: list of papers IDs (separated by ";") that reports the effect in the subtype for the marker

### Matrix
Options: ```-m FILENAME```/```--matrix-output FILENAME```
```
mutfinder -m mutations_output.tsv your_fasta.fa
mutfinder --matrix-output mutations_output.tsv your_fasta.fa
```

Outputs amino acids present for each sample in each position of any mutation found in the analysis.
The output file is a tab-separated matrix with samples on rows and mutations on columns.

## FASTA header parsing
To assign a FASTA sequence to a specific segment MutFinder relies on the FASTA header.
In the header must be present both a sample ID (consistent among segments of the same sample),
and the segment name as reported in the database (eg. for avian influenza possible segments are ```PB2```, ```PB1```, ```PA```, ```HA```, ```NP```, ```NA```, ```MP```, ```NS```).
By default the program expects the sample ID first and then the segment name separated by an undrescore and nothing else (eg. ```my_sample_PB2```)

If your sequences have a different pattern in your name you can specify it with the option ```-n```/```--name-regex```.
It works with 2 capture groups, the first one for the sample and the second one for the segment.
Named groups (respectively ```sample``` and ```segment```) can be used to invert groups order.
Here you can find some examples for the regular expression:
| FASTA header                  | regular expression                    |
| ----------------------------- | ------------------------------------- |
| my_sample_PB2                 | ```(.+)_(.+)```                       |
| my_sample\|PB2                | ```(.+)\\|(.+)```                     |
| my_sample-PB2                 | ```(.+)-(.+)```                       |
| my_sample-PB2_something_else  | ```(.+)-(.+?)_```                     |
| PB2_my_sample                 | ```(?P<segment>.+?)_(?P<sample>.+)``` |


- Option ```--skip-unmatch-names```: skip sequence and continue analysis when a FASTA header does not match the regular expression provided. By default exit with error.

- Option ```--skip-unknown-segments```: skip sequence and continue analysis when the segment found in FASTA header is not present in the database. By default exits with error. 

> **_NOTE:_** to find the regular expression that better fits your data you can try it on [Regex101](https://regex101.com/) selecting ```Python``` flavor.
