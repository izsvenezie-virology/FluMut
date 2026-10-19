<p align="center">
    <picture>
      <source media="(prefers-color-scheme: dark)" srcset="docs/images/flumut_logo_w.png">
      <source media="(prefers-color-scheme: light)" srcset="docs/images/flumut_logo_b.png">
      <img alt="FluMut logo." src="docs/images/flumut_logo_b.png">
    </picture>
</p>

[![GitHub Release](https://img.shields.io/github/v/release/izsvenezie-virology/FluMut?label=FluMut)](https://github.com/izsvenezie-virology/FluMut/releases/latest/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/flumut/README.html)
[![install with pip](https://img.shields.io/badge/install%20with-pip-brightgreen.svg)](https://pypi.org/project/flumut/)

FluMut is an open-source tool designed to search for molecular markers with potential impact on the biological characteristics of Influenza A viruses, starting from complete or partial nucleotide genome sequences.
It works on all HA (H1–H16) and NA (N1–N9) subtypes and is available both as a command line tool and as a graphical interface (FluMutGUI).

For the complete documentation please visit [FluMut site](https://izsvenezie-virology.github.io/FluMut/).

## Installation

### FluMutGUI

FluMutGUI is the graphical interface of FluMut, easy to use and easy to install.
If you need the CLI version, check [FluMut installation](#flumut).
Download the executable for your OS and double click on it to launch FluMutGUI, no installation required:

- [Windows](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/FluMutGUI-Windows.exe)
- [MacOS](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/FluMutGUI-macOS)
- [Linux](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/FluMutGUI-Linux)

FluMutGUI can also be installed from [PyPI](https://pypi.org/project/flumut/) (`pip install "flumut[gui]"`) and launched with `flumutgui`.

### FluMut

FluMut is a CLI tool, if you prefer a graphical interface see [FluMutGUI](#flumutgui) installation.

#### Pip [![install with pip](https://img.shields.io/badge/install%20with-pip-brightgreen.svg)](https://pypi.org/project/flumut/)

FluMut is available on [PyPI](https://pypi.org/project/flumut/).
Before installing FluMut via Pip you need:

- [Python](https://www.python.org/downloads/) 3.10 or newer
- [Pip](https://pypi.org/project/pip/) (often packed with Python)

Then, you can install FluMut with this simple command:

```
pip install flumut
```

#### Bioconda [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/flumut/README.html)

FluMut is also available on [Bioconda](https://bioconda.github.io/recipes/flumut/README.html).
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

### Input

FluMut can analyze multiple Influenza A virus sequences simultaneously.
It can handle partial and complete genome sequences of multiple samples.
You must provide one or more files containing the nucleotide sequences in FASTA format.
Sequences must adhere to the [IUPAC code](https://www.bioinformatics.org/sms/iupac.html).

FluMut relies on the FASTA header to assign the sequence to a sample.
For this reason, the header must contain a sample ID consistent among sequences of the same sample.
By default FluMut expects the sample ID followed by an underscore and the segment name (e.g. `my_sample_PB2`); a different format can be parsed with the `-n`/`--name-regex` option.
The segment is detected automatically by aligning the sequence against the references in the database: the one in the header is used only as a hint.
Find out more [here](https://izsvenezie-virology.github.io/FluMut/docs/usage/input-file).

An example of input file can be downloaded [here](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/fasta_input_example.fa).

#### Basic usage

You can get the output file in an [Excel format](#excel) (user-friendly) running:

```
flumut -x excel_output.xlsm your_fasta.fa
```

If you prefer the [text outputs](#text-outputs) (machine-readable format) run:

```
flumut -m markers_output.tsv -M mutations_output.tsv -l literature_output.tsv your_fasta.fa
```

You can analyze several FASTA files at once by listing them all at the end of the command.
The full list of options is available [here](https://izsvenezie-virology.github.io/FluMut/docs/usage/usage-cli#options).

#### Update database

The markers database is shipped inside FluMut, and every database update is released as a new FluMut version.
You can check whether a newer one is available with:

```
flumut --check-update
```

If you cannot update FluMut, you can download the latest database and use it with the `-D`/`--db-file` option:

```
flumut --download-db flumut_db.sqlite
flumut -D flumut_db.sqlite -x excel_output.xlsm your_fasta.fa
```

### Outputs

FluMut can produce an Excel output or text outputs:

- [Excel](#excel)
- [Text outputs](#text-outputs)

By default FluMut reports only markers where all mutations are found.
You can report all markers where at least one mutation is found using option `-r`/`--relaxed`.

#### Excel

This is the most user-friendly and complete output.
You can obtain this output using the `-x`/`--excel-output` option.
Find out more [here](https://izsvenezie-virology.github.io/FluMut/docs/output#excel-output).

The Excel output also contains a `Checks` sheet listing possible issues found in your sequences (frameshifts, premature stop codons, duplicated proteins, etc.).

> **_IMPORTANT:_** To enable the navigation feature the extension of the Excel file must be `.xlsm`.
> If you don't care about navigation, you can use `.xlsx` extension.
> Other extensions lead to unreadable files.

#### Text outputs

You can obtain 3 different text outputs:
| Option | Output | Description |
| -- | -- | -- |
| `-m`/`--markers-output` | [Markers output](https://izsvenezie-virology.github.io/FluMut/docs/output#markers-output) | List of detected markers |
| `-M`/`--mutations-output` | [Mutations output](https://izsvenezie-virology.github.io/FluMut/docs/output#mutations-output) | List of amino acids present in the positions of mutations of interest for each sample |
| `-l`/`--literature-output` | [Literature output](https://izsvenezie-virology.github.io/FluMut/docs/output#literature-output) | List of the papers describing the detected markers |

## Cite FluMut

If you use FluMut, please cite:

> Giussani, E., Sartori, A. et al. (2025). FluMut: a tool for mutation surveillance in highly pathogenic H5N1 genomes. Virus Evolution, [10.1093/ve/veaf011](https://doi.org/10.1093/ve/veaf011).

## License

FluMut is licensed under the GNU Affero v3 license (see [LICENSE](LICENSE)).

## Fundings

This work was partially supported by the FLU-SWITCH Era-Net ICRAD (grant agreement No. 862605), by EU funding under the NextGeneration EU-MUR PNRR Extended Partnership initiative on Emerging Infectious Diseases (Project No. PE00000007, INF-ACT), and by KAPPA-FLU HORIZON-CL6-2022-FARM2FORK-02-03 (grant agreement No. 101084171).

<p align="center" margin="10px">
    <img style="height:80px;margin:8px" alt="Logo supporter, FLU-SWITCH" src="docs/images/Logo-Flu-Switch.png"/>
    <img style="height:80px;margin:8px" alt="Logo supporter, INF-ACT" src="docs/images/Logo-Inf-act.jpg"/>
    <img style="height:80px;margin:8px" alt="Logo supporter, European Union" src="docs/images/Logo-eu.png"/>
    <img style="height:80px;margin:8px" alt="Logo supporter, KAPPA-FLU" src="docs/images/logo-kappa-flu.jpg"/>
</p>

> Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Health and Digital Executive Agency (HEDEA).
> Neither the European Union nor the granting authority can be held responsible for them
