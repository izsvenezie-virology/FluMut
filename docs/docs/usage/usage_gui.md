---
title: GUI usage
layout: default
parent: Usage
nav_order: 2
permalink: docs/usage/usage-gui
---

# Basic Usage
FluMutGUI is very simple to use:
1. Update the database to latest version
1. Select the FASTA file you want to analyze (learn more [here](./input-file))
1. Select which [outputs](../output) you want (and change the file name if you want)
1. Launch the analysis

![](../../images/GUI-usage.png)

FluMut will analyze your samples and will create the selected outputs.
When it finishes check the messages and then you can close the program.

![](../../images/GUI-usage-done.png)

{: .important}
>FluMutGUI has the [`--skip-unmatch-names`](./usage-cli#options) and [`--skip-unknown-segments`](./usage-cli#options) options flagged by default.
>
>Read the log for a list of skipped sequences.
