# Sketchy <a href='https://github.com/esteinig'><img src='img/logo.png' align="right" height="210" /></a>

![](https://img.shields.io/badge/version-alpha-red.svg)
![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![](https://img.shields.io/badge/docs-github-green.svg)
![](https://img.shields.io/badge/BioRxiv-v1-green.svg)

Real-time lineage matching and genotyping from uncorrected nanopore reads

## Overview

**`v0.3-alpha: test build, still no tests`**

`Sketchy` is an online lineage matching algorithm for real-time genotyping and susceptibility prediction in bacterial pathogens using nanopore sequencing platforms. `Sketchy` uses MinHash distances to match reads to a representative population genomic sketch of the target organism. Because it relies on well characterized, high-quality genomes from the public domain, it is currently restricted to common pathogens for which we have sufficient data that can be processed with the [`pf-core/pf-survey`](https://github.com/pf-core) Nextflow pipeline, which at this stage supports *Staphylococcus aureus*, *Klebsiella pneumoniae* and *Mycobacterium tuberculosis* specific typing. *Burkholderia pseudomallei* with antibiotic resistance genotype using `ArDAP` and other candidates are planned for beta release.

Preprint coming soon.

## Install

```
git clone https://github.com/esteinig/sketchy
pip install sketchy
```

