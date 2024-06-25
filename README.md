# FluMut

FluMut is an open-source tool designed to search for molecular markers with potential impact on the biological characteristics of Influenza A viruses of the A(H5N1) subtype, starting from complete or partial nucleotide genome sequences.

For the complete documentation please visit [FluMut site](https://izsvenezie-virology.github.io/FluMut/).

## Installation

### Prerequisites
FluMut is available for Windows, Linux and macOS.

#### Pip
FluMut is available on [PyPI](https://pypi.org/flumut).
Before installing FluMut via Pip you need:
- [Python](https://www.python.org/downloads/)
- [Pip](https://pypi.org/project/pip/) (often packed with Python)

Then, you can install FluMut with this command:
```
pip install flumut
```

#### Bioconda
FluMut is also available on [Bioconda](https://bioconda.github.io/flumut).
You can install using Conda or Mamba.
- [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (recommended)
```
mamba install -c bioconda flumut
```
- [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
```
conda install -c bioconda flumut
```

## Usage
### Basic usage
You can get the output file in an [Excel format](#excel) (user-friendly) running:
```
flumut -x excel_output.xlsm your_fasta.fa
```
If you prefer the [text outputs](#text-outputs) (machine-readable format) run:
```
flumut -m markers_output.tsv -M mutations_output.tsv -l literature_output.tsv your_fasta.fa
```

### Update database
You should always use the latest version of our database and you can do it just by running this command:
```
flumut --update
```

## Input
FluMut can analyze multiple A(H5N1) Influenza virus sequences simultaneously.
It can handle partial and complete genome sequences of multiple samples.
You must provide a single file containing all the nucleotide sequences in FASTA format.
Sequences must adhere to the [IUPAC code](https://www.bioinformatics.org/sms/iupac.html).

FluMut relies on the FASTA header to assign the sequence to a specific segment and sample.
For this reason, the header must contain both a sample ID (consistent among sequences of the same sample) and one of the the following segment names: `PB2`, `PB1`, `PA`, `HA`, `NP`, `NA`, `MP`, `NS`.

An example of input file can be downloaded [here](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/fasta_input_example.fa).

## Outputs
FluMut can produce an Excel output or text outputs:
- [Excel](#excel)
- [Text outputs](#text-outputs)

By default FluMut reports only markers where all mutations are found.
You can report all markers where at least one mutation is found using option `-r`/`--relaxed`.

### Excel
This is the most user-friendly and complete output. 
You can obtain this output using the `-x`/`--excel-output` option.
Find out more [here](https://izsvenezie-virology.github.io/FluMut/docs/output#excel-output).

>**_IMPORTANT:_** To enable the navigation feature the exstension of the Excel file must be `.xlsm`.
>If you don't care about navigation, you can use `.xlsx` exstension.
>Other exstensions lead to unreadable files.

### Text outputs
You can obtain 3 different text outputs:
| Option | Output | Desctription |
| -- | -- | -- |
| `-m`/`--markers-output` | [Markers output](https://izsvenezie-virology.github.io/FluMut/docs/output#markers-output) | List of detected markers |
| `-M`/`--mutations-output` | [Mutations output](https://izsvenezie-virology.github.io/FluMut/docs/output#mutations-output) | List of amino acids present in the positions of mutations of interest for each sample |
| `-l`/`--literature-output` | [Literature output](https://izsvenezie-virology.github.io/FluMut/docs/output#literature-output) | List of all papers present in the database |

## Cite FluMut
We are currently writing the paper. 
Until the publication please cite the GitHub repository:

[https://github.com/izsvenezie-virology/FluMut](https://github.com/izsvenezie-virology/FluMut)

## License
FluMut is licensed under the GNU Affero v3 license (see [LICENSE](LICENSE)).

# Fundings

This work was supported by FLU-SWITCH Era-Net ICRAD (grant agreement No 862605) and by the NextGeneration EU-MUR PNRR Extended Partnership initiative on Emerging Infectious Diseases (Project no. PE00000007, INF-ACT)

![](docs/images/Logo-Flu-Switch.png) ![](docs/images/Logo-Inf-act.jpg) ![](docs/images/Logo-eu.png)

>Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Health and Digital Executive Agency (HEDEA). 
>Neither the European Union nor the granting authority can be held responsible for them
