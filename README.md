# sketchy <a href='https://github.com/esteinig'><img src='img/logo.png' align="right" height="210" /></a>


![](https://img.shields.io/badge/lang-rust-black.svg)
![](https://img.shields.io/badge/version-0.4.0-purple.svg)
![](https://img.shields.io/badge/biorxiv-v1-blue.svg)

Real-time lineage matching and genotyping from uncorrected nanopore reads

## Overview

**`v0.4.0: beta, rust core libs`**

`Sketchy` is an online lineage calling and genotyping algorithm based on the heuristic principle of genomic neighbor typing by [Karel Břinda and colleagues (2019)](https://www.biorxiv.org/content/10.1101/403204v2). `Sketchy` computes the sum of min-wise hashes shared with species-wide sketches of bacterial pathogen genomes and their associated genotypes e.g. multi-locus sequence types, susceptibility profiles computed with [Mykrobe](https://github.com/Mykrobe-tools/mykrobe) or serotype alleles inferred with [Kleborate](https://github.com/katholt/kleborate). A list of precomputed genotype features can be found in the corresponding pathogen reference  sections.

Currently supported species are:

* *Staphylococcus aureus*
* *Klebsiella pneumoniae* 
* *Mycobacterium tuberculosis*

Please see our preprint for guidance on the limitations of `Sketchy`.

- [Install](#install)
  - [:snake: `conda`](#conda)
  - [:whale: `docker`](#docker)
  - [:new_moon: `singularity`](#singularity)
  - [:rocket: `cargo`](#cargo)
- [Setup](#setup)
- [Usage](#usage)
  - [Python CLI](#python-client)
  - [Rust CLI](#rust-client)
  - [Evaluation outputs](#rust-client)
  - [Online sequencing run](#rust-client)
  - [Android mobile phones](#rust-client)
- [How it works](#how-it-works)
- [Reference sketches](#reference-sketches)
  - [*Staphylococcus aureus*](#rust-client-tasks)
  - [*Klebsiella pneumoniae*](#rust-client-tasks)
  - [*Mycobacterium tuberculosis*](#rust-client-tasks)
- [Constructing reference sketches](#reference-sketches)
  - [Genome assemblies and sketch construction](#rust-client-tasks)
  - [Genotype features and index preparation](#rust-client-tasks)
  - [Lineage and local sketches](#rust-client-tasks)
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

Sketchy implements a `Rust` command-line interface (`sketchy-rs`) for computation and evaluation on read streams and a `Python` command-line interface (`sketchy`) for evaluation plots and other utilities. It is recommended to use one of the following options to install the required dependencies and access the complete computation and evaluation pipeline.

#### `Conda`

`Sketchy` is currently on a private channel and requires some dependencies from `BioConda`:

```sh
conda install -c esteinig -c bioconda sketchy
```

#### `Docker`

The `Docker` container is based on the `Alpine` image with internal `Conda` environments:

```sh
docker pull esteinig/sketchy
docker run esteinig/sketchy --help
```

#### `Singularity`

`Docker` containers work with `Singularity`:

```sh
singularity exec docker://esteinig/sketchy sketchy --help
```

#### `Cargo`

For the pure `Rust` client, where the `compute` subtask requires `Mash` in `$PATH`:

```bash
cargo install sketchy-rs
```

On `Linux` systems one should be able to install `Mash` conveniently e.g.

```bash
sudo apt-get install mash
```


## Setup

Pull default species sketches into local storage before first use:

```
sketchy database pull
```

Local sketches and template names can be viewed with:

```
sketchy database list
```

## Usage

See the `Tasks and Parameters` section for details on all tasks and settings available in `Sketchy`.

### Python CLI

`Sketchy` can be run through a wrapper in the Python CLI which is suitable for completed read files. Read streams and online sequencing runs should be served with the `Rust` CLI.

```bash
sketchy run --help
```

Species-wide sketch templates are available for:

* *Staphylococcus aureus*: `saureus`
* *Klebsiella pneumoniae*: `kpneumoniae`
* *Mycobacterium tuberculosis*: `mtuberculosis`

When using a template, execution looks like this:

```bash
sketchy run \
  --fastx     test.fq \
  --reads     5000 \
  --ranks     20 \
  --sketch    saureus \
  --prefix    path/to/prefix \
  --stable    1000 \
  --palette   YnGnBu
```

**Custom sketches**

Custom reference sketch collections can be generated with the [`sketchy feature prepare`](#sketchy-feature-prepare) task as described in the [`Constructing reference sketches`](#constructing-reference-sketches) section, and must include a:

* reference sketch
* numeric genotype index
* genotype key file 

with the same file names and the following extensions, such that:

```bash
ref.msh   # sketch
ref.tsv   # index
ref.json  # key 
```

Custom collections can be used with reference to the path and file name in the `--sketch` option:

```bash
sketchy run \
  --fastx   test.fq \
  --reads   5000 \
  --sketch  path/to/ref
```

### Rust CLI

The `Rust` command line interface implements two subtasks: `sketchy-rs compute` and `sketchy-rs evaluate`. Both read from `/dev/stdin` and can be piped. Setup with `sketchy database pull` deposited the default sketches to `~/.sketchy/sketches` so we can set an environmental variable for convenience:
 
```bash
SKPATH=~/.sketchy/sketches
```
 
`Compute` internally calls `Mash` and processes the output stream by computing the sum of shared hashes. If heatmaps should be included in the evaluations, the output should be directed to a file, e.g.
 
 ```bash
 cat test.fq | head -20000 | \
 sketchy-rs compute \
    --sketch $SKPATH/saureus.msh \
    --ranks 20 \
 > test.ssh.tsv
 ```
 
`Evaluate` then computes the sum of ranked sums of shared hashes, and other summaries for plotting:

```bash
cat test.ssh.tsv | \
sketchy-rs evaluate \
    --index $SKPATH/saureus.tsv \
    --stable 1000 \
> test.sssh.tsv
```

The `Rust` pipeline can therefore be executed as:

```bash
cat test.fq | head -20000 \
| sketchy-rs compute \
    --sketch $SKPATH/saureus.msh \
    --ranks 20 \
| sketchy-rs evaluate \
    --index $SKPATH/saureus.tsv \
    --stable 1000 \
> test.sssh.tsv
```

Plotting and evaluation summaries are handled in the `Python` CLI and accessed via the `sketchy plot` task:

```
sketchy plot \
    --sssh test.sssh.tsv \
    --ssh test.ssh.tsv \
    --index $SKPATH/saureus.tsv \
    --key $SKPATH/saureus.json \
    --stable 1000 \
    --palette YnGnBu \
    --prefix test \
    --format png
```
 
### Sketchy evaluation outputs

### Online sequencing run

### Android mobile phones

To set up the `Rust` CLI on Android mobile phones, the following can be done in a couple of minutes:

1. Install the [`UserLAnd`](https://github.com/CypherpunkArmory/UserLAnd) app 
2. Setup an `Ubuntu` image
3. Run the following script

```bash
sudo apt-get update
sudo apt-get install curl mash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install sketchy-rs
```

Python CLI has not been tested.

## How it works

## Constructing reference sketches

Reference sketches can be rapidly constructed and prepared for use with `Sketchy`. Custom sketches are useful for prediction on species currently not offered in the default collection, lineage specific sub-sketches of a species, or local genome collections, such as from  healthcare providers or surveillance programs that are not publicly available. All that is needed is a set of high-quality reference assemblies and their associated genotype features. Ultimately, genome and feature representation are the most important components to consider, as they define the genomic neighbors that can be typed with `Sketchy`. 

### Genome assemblies and sketch construction

Assemblies should be of sufficient quality for genotyping and can produced e.g. with tools from the [`Torstyverse`](https://github.com/tseemann) like [`Shovill`](https://github.com/tseemann/shovill) or as part of large-scale public archive surveillance pipelines like [`Pathfinder`](https://github.com/pf-core). 

Given a set of high-quality assemblies in the current directory:

```
DRR083589.fasta
DRR083590.fasta
DRR119226.fasta
DRR119227.fasta
DRR128207.fasta
DRR128208.fasta
```

`Mash` can be used directly to construct the sketch with the following default parameters:

```
mash sketch -s 1000 -k 15 *.fasta
```

When constructing sketches from thousands of genomes, it might be more convenient to use the `Nextflow` pipeline, which parallelizes the sketch construction using `mash sketch` executions and `mash paste`:

```
nextflow run esteinig/sketchy --build true --fasta "*.fasta" --sketch_size 1000 --kmer_size 15 
```

See the [`Nextflow`](#nextflow) section for additional setting and the [`Benchmarks`](#benchmarks) section for guidance on selecting an appropriate sketch and k-mer size for `Sketchy`. 

### Genotype features and index preparation
