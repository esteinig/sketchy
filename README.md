# sketchy <a href='https://github.com/esteinig'><img src='docs/imalogo.png' align="right" height="250" /></a>

![](https://img.shields.io/badge/lang-rust-black.svg)
![](https://img.shields.io/badge/version-0.6.0-purple.svg)
![](https://img.shields.io/badge/biorxiv-1.0-blue.svg)

Genomic neighbour typing for lineage and genotype inference of bacterial pathogens

## Overview

**`v0.6.0`**

`Sketchy` is a lineage calling and genotyping tool based on the heuristic principle of genomic neighbor typing developed by [Karel Břinda and colleagues (2020)](https://www.biorxiv.org/content/10.1101/403204v2). It queries species-wide ('hypothesis-agnostic') reference sketches using MinHash and infers associated genotypes based on the closest match, including multi-locus sequence types, susceptibility profiles, virulence factors or other genome-associated features provided by the user. Unlike the original implementation, `Sketchy` does not use phylogenetic trees (which has some downsides) and is easier to modify for users, for example to construct reference sketches for local genome collections, or species not implemented in the current release.

See the [latest docs](https://esteinig.github.io/sketchy) for install, usage and database building.


## Preprints

Preprint assessing `Sketchy` for bacterial genotype predictions from low read numbers, using simulated data and novel nanopore data from an *S. aureus* outbreak in remote communities of Far North Queensland and Papua New Guinea. Metagenomic application in a cystic fibrosis patient undergoing antimicrobial therapy uses data from Tania Duarte et al. (2022) metagenomic surveillance study.

* Sketchy preprint:
* FNQ and PNG bioproject:
* Cystic fibrosis preprint:


## Citation

If you use `Sketchy` for research and other applications, please cite:

```

```
