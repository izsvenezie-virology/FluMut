---
title: Installation
layout: default
nav_order: 2
permalink: docs/installation
---

# Installation

FluMut comes in two flavours that share the same analysis engine and the same database:

- [FluMutGUI](#flumutgui): the graphical interface, easy to use and easy to install
- [FluMut](#flumut): the command line tool, best suited for pipelines and automation

## FluMutGUI

FluMutGUI is the graphical interface of FluMut.
If you need the CLI version, check [FluMut installation](#flumut).

### Standalone executable

Download the executable for your OS from the [latest FluMut release](https://github.com/izsvenezie-virology/FluMut/releases/latest/) and double click on it to launch FluMutGUI.
No installation is required.

- [Windows](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/FluMutGUI-Windows.exe)
- [MacOS](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/FluMutGUI-macOS)
- [Linux](https://github.com/izsvenezie-virology/FluMut/releases/latest/download/FluMutGUI-Linux)

{: .note}
On Linux and MacOS you may need to make the downloaded file executable (`chmod +x FluMutGUI-Linux`) before launching it.

### With Pip

FluMutGUI is also included in the FluMut Python package as an optional extra.
Follow the [Pip installation](#pip) instructions installing `flumut[gui]` instead of `flumut`:

```
pip install "flumut[gui]"
```

Then launch the graphical interface with:

```
flumutgui
```

## FluMut

FluMut is a CLI tool, if you prefer a graphical interface see [FluMutGUI](#flumutgui) installation.

### Pip [![install with pip](https://img.shields.io/badge/install%20with-pip-brightgreen.svg)](https://pypi.org/project/flumut/)

FluMut is available on [PyPI](https://pypi.org/project/flumut/).
Before installing FluMut via Pip you need:

- [Python](https://www.python.org/downloads/) 3.10 or newer
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

## Updates

The markers database is shipped inside FluMut, so a new database is distributed as a new FluMut release.
To get the latest markers, keep FluMut up to date:

- FluMutGUI shows a notice in the main window when a newer version is available
- with the CLI, run `flumut --check-update`

Then update FluMut the same way you installed it (e.g. `pip install --upgrade flumut` or download the new executable).

If you cannot update FluMut, you can still download the latest database and use it with your current installation, as described in [CLI usage](usage/usage-cli#update-the-database).
