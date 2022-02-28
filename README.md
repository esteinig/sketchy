# sketchy <a href='https://github.com/esteinig'><img src='docs/images/logo.png' align="right" height="250" /></a>

![](https://img.shields.io/badge/lang-rust-black.svg)
![](https://img.shields.io/badge/version-0.6.0-black.svg)
![](https://img.shields.io/badge/biorxiv-0.11.2-green.svg)

Genomic neighbor typing for lineage and genotype inference

## Overview

**`v0.6.0`**

`Sketchy` is a lineage calling and genotyping tool based on the heuristic principle of genomic neighbor typing developed by [Karel Břinda and colleagues (2020)](https://www.biorxiv.org/content/10.1101/403204v2). It queries species-wide ('hypothesis-agnostic') reference sketches using MinHash and infers associated genotypes based on the closest match, including multi-locus sequence types, susceptibility profiles, virulence factors or other genome-associated features provided by the user. Unlike the original implementation in [`RASE`](https://github.com/c2-d2/rase-pipeline), `Sketchy` does not use phylogenetic trees which has some downsides, e.g. for sublineage genotype predictions (see below). 

See the [latest docs](https://esteinig.github.io/sketchy) for install, usage and database building.

## Strengths and limitations

Please see the preprint for detailed limitations of `Sketchy`. 

* `Sketchy` performs best on lineage predictions and lineage-wide genotypes from few reads (you can try with < 1000 reads)
* reference sketches and genotype indices can be constructed easily from large genome collections (e.g. species sketches)

However:

* clade-specific genotype resolution is not as good as when using phylogenetic guide trees (e.g. with `RASE`) - for sublineage genotyping i.e. when the number of genomes in the reference sketch is smaller, sketch size (`s`) can be increased to compensate for some of the approximate matching (~ 10k - 100k)

> ⚠️ We are currently working on resolving this with reference sketch k-mer size optimisation

* `Sketchy` genotype inference may be difficult for species with high rates of homologous recombination - it is advised to run simulations and validations for species we have not provided reference sketches for

## Release v0.6.0 

Version used in preprint available from `master` branch or release `0.5.9` (implements all functions of `0.6.0` necessary for preprint):

```
git clone https://github.com/esteinig/sketchy
cd sketchy && cargo build --release
```

Target: 7 March

* Bioconda recipe update
* Nextflow pipelines for optimisation and validation
* Upload of species reference sketches
* Utility tools to build sketches with data from [Grace Blackwell et al.](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001421)
* Default k-mer sizes for the databases are currently very small (k = 16) - optimised databases will be available

## Preprint

If you use `Sketchy` for research and other applications, please cite:

>  Steinig et al. (2022) - Genomic neighbor typing for bacterial outbreak surveillance - bioRxiv 2022.02.05.479210; doi: https://doi.org/10.1101/2022.02.05.479210 
