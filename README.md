# sketchy <a href='https://github.com/esteinig'><img src='img/logo.png' align="right" height="210" /></a>

![](https://img.shields.io/badge/version-alpha-red.svg)
![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![](https://img.shields.io/badge/core-rust-black.svg)
![](https://img.shields.io/badge/docs-github-green.svg)
![](https://img.shields.io/badge/BioRxiv-v1-orange.svg)

Real-time lineage matching and genotyping from uncorrected nanopore reads

## Overview

**`v0.4.0: public build, rust core`**

`Sketchy` is an online lineage calling and genotyping algorithm for bacterial pathogens using uncorrected nanopore reads. Currently supported species are *Staphylococcus aureus*,  *Klebsiella pneumoniae* and *Mycobacterium tuberculosis*. Please see the preprint for guidance on the limitations of `sketchy`.

- [Install](#install)
  - [`conda`](#conda)
  - [`docker`](#docker)
  - [`singularity`](#singularity)
  - [`cargo`](#cargo)
- [Setup](#setup)
- [Usage](#usage)
  - [Basic usage](#basic-usage)
- [How it works](#how-it-works)
- [Citing](#citing)
  - [Bibtex](#bibtex)

## Install

Sketchy implements a `Rust` command-line interface (`sketchy-rs`) for computation on read streams and a `Python` command-line interface (`sketchy`) for evaluation plots and other utilities. It is recommended to use one of the following options to install the required dependencies and access the complete computation and evaluation pipeline.

#### `Conda`

`Sketchy` is currently on a private channel and requires some dependencies from `BioConda`:

```sh
conda install -c esteinig -c bioconda sketchy
```

#### `Docker`

The Docker container is based on the `Alpine` image with internal `Conda` environments:

```sh
docker pull esteinig/sketchy
docker run esteinig/sketchy sketchy --help
```

For a lean `Alpine` container containing only the `Rust` implementation of `sketchy` use:

```sh
docker pull esteinig/sketchy-rs
docker run esteinig/sketchy-rs -h
```

#### `Singularity`

You can use the `Docker` containers with `Singularity`:

```sh
singularity exec docker://esteinig/sketchy sketchy --help
```

#### `Cargo`

For the bare-bones Rust libraries without evaluation plots:

```sh
cargo install sketchy-rs
```

## Setup

Pull default species sketches into local storage:

`sketchy database pull`

Local sketches can be viewed with:

`sketchy database list`

## Usage

### Basic usage

### How it works
