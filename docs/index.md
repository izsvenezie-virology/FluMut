---
title: Home
layout: home
nav_order: 1
permalink: /
---

![](images/flumut_logo_b.svg)

[![GitHub Release](https://img.shields.io/github/v/release/izsvenezie-virology/FluMut?label=FluMut)](https://github.com/izsvenezie-virology/FluMut/releases/latest/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/flumut/README.html)
[![install with pip](https://img.shields.io/badge/install%20with-pip-brightgreen.svg)](https://pypi.org/project/flumut/)

FluMut is an open-source tool designed to search for molecular markers with potential impact on the biological characteristics of avian Influenza A viruses, starting from complete or partial nucleotide genome sequences.

A marker is defined as a set of amino acid mutations (one or more mutations located in the same or different proteins) for which a phenotypic effect has been demonstrated.
FluMut is based on a database containing a list of markers previously described in the literature as associated with:

- increased zoonotic potential
- increased virulence in different host categories
- increased/conferred resistance towards antiviral drugs

This database is periodically updated with newly described markers, and each FluMut release ships the latest version of it.

For each marker identified within the submitted sequences, the tool provides information on:

- the position of the mutations composing the marker
- the effect of the marker as described in the literature
- the influenza subtype in which the marker was observed or tested
- the corresponding bibliographic references

FluMut aligns each sequence against a set of reference sequences, so it works on all Influenza A subtypes (HA H1–H16 and NA N1–N9) and can recognize the segment of each sequence by itself.
While translating the sequences, it also reports possible quality issues such as frameshifts, premature stop codons and duplicated segments.

In addition, the tool generates an Excel workbook for a better visualisation of the results, as described in the [outputs section](docs/output#excel-output).

FluMut is available both as a command line tool and as a graphical interface (FluMutGUI); see the [installation page](docs/installation).
