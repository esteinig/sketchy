# sketchy <a href='https://github.com/esteinig'><img src='docs/images/logo.png' align="right" height="250" /></a>

![](https://img.shields.io/badge/lang-rust-black.svg)
![](https://img.shields.io/badge/version-0.6.0-purple.svg)
![](https://img.shields.io/badge/biorxiv-1.0-blue.svg)

Genomic neighbour typing for lineage and genotype inference of bacterial pathogens

## Overview

**`v0.6.0`**

`Sketchy` is a lineage calling and genotyping tool based on the heuristic principle of genomic neighbor typing developed by [Karel Břinda and colleagues (2020)](https://www.biorxiv.org/content/10.1101/403204v2). It queries species-wide ('hypothesis-agnostic') reference sketches using MinHash and infers associated genotypes based on the closest match, including multi-locus sequence types, susceptibility profiles, virulence factors or other genome-associated features provided by the user. Unlike the original implementation in [`RASE`](https://github.com/c2-d2/rase-pipeline), `Sketchy` does not use phylogenetic trees which has some downsides, e.g. for sublineage genotype predictions (see below). 

See the [latest docs](https://esteinig.github.io/sketchy) for install, usage and database building.

## Strengths and limitations

Please see the preprint for detailed limitations of `Sketchy`. 

* `Sketchy` performs best on lineage predictions and lineage-distributed genotypes from few reads
* Reference sketches and genotype indices can be constructed easily from large genome collections

However:

* Because of the approximate matching approach using MinHash, sub-lineage genotype resolution is not as good as `RASE`, which uses phylogenetic trees (and may be preferred for inference of clade-specific traits)
* `Sketchy` genotype inference may be difficult for species with high rates of homologous recombination - it is advised to run simulations and validations for species we have not provided reference sketches

## Preprints

Preprint assessing `Sketchy` for bacterial genotype predictions from low read numbers, using simulated data and novel nanopore data from an *S. aureus* outbreak in remote communities of Far North Queensland and Papua New Guinea. 

* Sketchy preprint:
* FNQ and PNG preprints:


## Citation

If you use `Sketchy` for research and other applications, please cite:

```

```
