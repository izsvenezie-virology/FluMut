---
title: Installation
layout: default
nav_order: 2
permalink: docs/installation
---

# Installation

## FluMut-GUI

FluMut-GUI is the graphical interface of FluMut, easy to use and easy to install.
If you need the CLI version, check [FluMut installation](#flumut).
Download the compiled version for your OS, double click on it and follow instructions to install FluMutGUI:

- [Windows](https://github.com/izsvenezie-virology/FluMutGUI/releases/latest/download/FluMutGUI-Windows-installer.exe)
- [MacOS](https://github.com/izsvenezie-virology/FluMutGUI/releases/latest/download/FluMutGUI-MacOS.zip)
- [Linux](https://github.com/izsvenezie-virology/FluMutGUI/releases/latest/download/FluMutGUI-linux.tar.gz)

## FluMut

FluMut is a CLI tool, if you prefer a graphical interface see [FluMut-GUI](#flumut-gui) installation.

### Pip [![install with pip](https://img.shields.io/badge/install%20with-pip-brightgreen.svg)](https://pypi.org/project/flumut/)

FluMut is available on [PyPI](https://pypi.org/project/flumut/).
Before installing FluMut via Pip you need:

- [Python](https://www.python.org/downloads/)
- [Pip](https://pypi.org/project/pip/) (often packed with Python)

Then, you can install FluMut with this simple command:

```
pip install flumut
```

{: .important}
When installing via Pip it's strongly recommended to use a [virtual environment](https://docs.python.org/3/library/venv.html).

### Bioconda [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/flumut/README.html)

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
