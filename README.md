# sketchy <a href='https://github.com/esteinig'><img src='docs/logo.png' align="right" height="250" /></a>

![](https://img.shields.io/badge/lang-rust-black.svg)
![](https://img.shields.io/badge/version-0.6.0-purple.svg)
![](https://img.shields.io/badge/biorxiv-1.0-blue.svg)

Genomic neighbour typing for lineage and genotype inference of bacterial pathogens

## Overview

**`v0.6.0`**

`Sketchy` is a lineage calling and genotyping tool based on the heuristic principle of genomic neighbor typing developed by [Karel Břinda and colleagues (2020)](https://www.biorxiv.org/content/10.1101/403204v2). It queries species-wide ('hypothesis-agnostic') reference sketches using MinHash and infers associated genotypes based on the closest match, including multi-locus sequence types, susceptibility profiles, virulence factors or other genome-associated features provided by the user. Unlike the original implementation, `Sketchy` does not use phylogenetic trees (which has some downsides) and is easier to modify for users, for example to construct reference sketches for local genome collections, or species not implemented in the current release.

- [Install](#install)
- [Usage](#usage)
- [Citing](#citing)

## Install

Sketchy implements a Rust command-line interface (`sketchy`) for computation and evaluation on read streams and a Python command-line application (`sketchy-utils`) for database construction and validation.

### `Cargo`

```bash
cargo install sketchy-rs
```

### `BioConda`

`Sketchy` is available on `BioConda` (thanks to [@mbhall88](https://github.com/mbhall88))

```
conda install -c conda-forge -c bioconda sketchy=0.5.0
```

### `Binaries`

Binaries for Linux and MacOS are available from the latest release.

