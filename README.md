# sketchy <a href='https://github.com/esteinig'><img src='img/logo.png' align="right" height="210" /></a>


![](https://img.shields.io/badge/lang-rust-black.svg)
![](https://img.shields.io/badge/version-0.4.0-purple.svg)
![](https://img.shields.io/badge/biorxiv-v1-blue.svg)

Real-time lineage matching and genotyping from uncorrected nanopore reads

## Overview

**`v0.4.0: beta, rust core libs`**

`Sketchy` is an online lineage calling and genotyping algorithm based on the heuristic principle of genomic neighbor typing by [Karel Břinda and colleagues (2019)](https://www.biorxiv.org/content/10.1101/403204v2). `Sketchy` computes the sum of min-wise hashes shared with species-wide sketches of bacterial pathogen genomes and their associated genotypes, for example multi-locus sequence types, susceptibility profiles computed with [Mykrobe](https://github.com/Mykrobe-tools/mykrobe) or serotype alleles inferred with [Kleborate](https://github.com/katholt/kleborate).

Currently supported species are:

* *Staphylococcus aureus*
* *Klebsiella pneumoniae* 
* *Mycobacterium tuberculosis*

Please see our preprint for guidance on the limitations of `Sketchy`.

- [Install](#install)
  - [`conda`](#conda)
  - [`docker`](#docker)
  - [`singularity`](#singularity)
  - [`cargo`](#cargo)
- [Setup](#setup)
- [Usage](#usage)
  - [Python CLI](#python-client)
  - [Rust CLI](#rust-client)
- [How it works](#how-it-works)
- [Tasks and parameters](#tasks)
  - [Rust CLI](#rust-client-tasks)
    - [`sketchy-rs compute`](#sketchy-rust-compute)
    - [`sketchy-rs evaluate`](#sketchy-rust-evaluate)
  - [Python CLI](#python-client-tasks)
    - [`sketchy run`](#sketchy-run)
    - [`sketchy plot`](#sketchy-plot)
    - [`sketchy utils`](#sketchy-utils)
    - [`sketchy feature`](#sketchy-feature)
    - [`sketchy database`](#sketchy-database)
- [Nextflow subworkflows](#nextflow)
  - [`parallel sketchy`](#nextflow-sketchy)
  - [`parallel sketch`](#nextflow-sketch)
  - [`metagenome`](#nextflow-metagenome)
- [Citing](#citing)
  - [BioRxiv](#bioarxiv)
  - [BibTeX](#bibtex)

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

For a lean `Rust` implementation of `Sketchy` based on `Alpine` use:

```sh
docker pull esteinig/sketchy-rs
docker run esteinig/sketchy-rs sketchy compute -h
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
