# Sketchy <a href='https://github.com/esteinig'><img src='img/logo.png' align="right" height="250" /></a>

Real-time lineage matching and genotyping from uncorrected nanopore reads

![](https://img.shields.io/badge/version-alpha-red.svg)
![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![](https://img.shields.io/badge/docs-latest-green.svg)
![](https://img.shields.io/badge/BioRxiv-incoming-green.svg)

## Overview

**`v0.1-alpha: internal pre-release, no tests`**

`Sketchy` is an online lineage matching algorithm for real-time typing of common pathogens using nanopore technology. It is a variant of inexact lineage matching as implemented by [Brinda et al. (2019)](https://www.biorxiv.org/content/early/2018/08/29/403204). `Sketchy` uses MinHash distances to match reads to a representative population genomic sketch of the target organism. Because it relies on well characterized, high-quality genomes from the public domain, it is currently restricted to common pathogens for which we have sufficient data, such as *Staphylococcus aureus* or *Klebsiella pneumoniae*.

Preprint coming soon.



