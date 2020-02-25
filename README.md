# sketchy <a href='https://github.com/esteinig'><img src='docs/logo.png' align="right" height="210" /></a>


![](https://img.shields.io/badge/lang-rust-black.svg)
![](https://img.shields.io/badge/version-0.4.2-purple.svg)
![](https://img.shields.io/badge/biorxiv-v1-blue.svg)

Real-time lineage matching and genotyping from uncorrected nanopore reads

## Overview

**`v0.4.2: beta, rust core libs`**

`Sketchy` is an online lineage calling and genotyping algorithm based on the heuristic principle of genomic neighbor typing by [Karel Břinda and colleagues (2020)](https://www.biorxiv.org/content/10.1101/403204v2). `Sketchy` computes the sum of min-wise hashes shared with species-wide reference sketches of bacterial pathogen genomes and their associated genotypes e.g. multi-locus sequence types, susceptibility profiles computed with [Mykrobe](https://github.com/Mykrobe-tools/mykrobe) or serotype alleles inferred with [Kleborate](https://github.com/katholt/kleborate). A list of precomputed genotype features can be found in the corresponding pathogen reference  sections.

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
  - [Python command line](#python-client)
  - [Rust command line](#rust-client)
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
  - [Python CLI](#python-client-tasks)
- [Nextflow pipeline](#nextflow)
- [Citing](#citing)
  - [BioRxiv](#bioarxiv)
  - [BibTeX](#bibtex)

## Install

Sketchy implements a `Rust` command-line interface (`sketchy-rs`) for computation and evaluation on read streams and a `Python` command-line interface (`sketchy`) for evaluation plots and other utilities. It is recommended to use one of the following options to install the required dependencies and access the complete computation and evaluation pipeline.


#### `Singularity`

I prefer `Singularity` for integrated access to host sytem files:

```sh
singularity pull docker://esteinig/sketchy:latest
./sketchy_latest.sif --help
```

#### `Docker`

`Docker` is ok too:

```sh
docker pull esteinig/sketchy:latest
docker run -it esteinig/sketchy --help
```

But to share files between container and the host system you need to set bindmounts, e.g. link the current directory with a hypothetical `test.fq` file to the preconfigured working directory `/data` inside the container using the current user permissions:

```sh
docker run -it \
  -v $(pwd):/data \
  -u $(id -u):$(id -g) \
  esteinig/sketchy run \
  --fastq /data/test.fq \
  --output /data/test
```


#### `Conda`

NOT AVAILABLE YET - `Sketchy` is currently on a private channel and requires some dependencies from `BioConda`:

```sh
conda install -c esteinig -c bioconda sketchy
```

#### `Cargo`

NOT AVAILABLE YET - `Rust` client only, where the `compute` subtask requires `Mash` in the system's `$PATH`:

```bash
cargo install sketchy-rs
```

On `Linux` systems one should be able to install `Mash` conveniently e.g.

```bash
sudo apt-get install mash
```


## Setup

Pull default species sketches into local storage before first use. ~ 2 GB for multiple species and sketch sizes:

```
sketchy pull
```

Local sketches and template names can be viewed with:

```
sketchy list
```

## Usage

See the `Tasks and Parameters` section for details on all tasks and settings available in `Sketchy`.

### Python CLI

`Sketchy` can be run through a wrapper in the `Python CLI` which is suitable for completed read files. Read streams and online sequencing runs should be served with the `Rust CLI`.

```bash
sketchy run --help
```

Species-wide sketch templates are available for:

* *Staphylococcus aureus*: `saureus`
* *Klebsiella pneumoniae*: `kpneumoniae`
* *Mycobacterium tuberculosis*: `mtuberculosis`

When using a template, execution looks like this:

```bash
sketchy run --fastq test.fq --sketch saureus --ranks 20 --limit 5000
```

More options can be viewed with

```bash
sketchy run --help
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
sketchy run --fastq test.fq --sketch ref
```

### Rust CLI

The `Rust` command line interface implements two subtasks: `sketchy-rs compute` and `sketchy-rs evaluate`. Both read from `/dev/stdin` and can be piped. Setup with `sketchy pull` deposited the default sketches to `~/.sketchy` so we can set an environmental variable for convenience:
 
```bash
SKPATH=~/.sketchy/saureus
```
 
`Compute` internally calls `Mash` and processes the output stream by computing the sum of shared hashes. If heatmaps should be included in the evaluations, the output should be directed to a file, e.g.
 
 ```bash
 cat test.fq | head -20000 | \
 sketchy-rs compute \
    --sketch $SKPATH/saureus.msh \
    --ranks 20 \
    --progress 1 \
    --threads 4 \
 > test.ssh.tsv
 ```
 
`Evaluate` then computes the sum of ranked sums of shared hashes, and other summaries for plotting:

```bash
cat test.ssh.tsv | \
sketchy-rs evaluate \
    --features $SKPATH/saureus.tsv \
    --stable 1000 \
> test.sssh.tsv
```

The `Rust` pipeline can be executed in one step, such as:

```bash
cat test.fq | head -20000 \
| sketchy-rs compute \
    --sketch $SKPATH/saureus.msh \
    --ranks 20 \
    --progress 1 \
    --threads 4 \
| sketchy-rs evaluate \
    --features $SKPATH/saureus.tsv \
    --stable 1000 \
> test.sssh.tsv
```

Plotting and evaluation summaries are handled in the `Python CLI` and accessed via the `sketchy plot` task:

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

To set up the `Rust CLI` on Android mobile phones, the following can be done in a couple of minutes:

1. Install the [`UserLAnd`](https://github.com/CypherpunkArmory/UserLAnd) app 
2. Setup an `Ubuntu` image
3. Run the following script

```bash
sudo apt-get update && sudo apt-get install curl mash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install sketchy-rs
wget https://storage.googleapis.com/np-core-sketchy/saureus.tar.gz \
  -O saureus.tar.gz && tar xvzf saureus.tar.gz
```

Reference sketch collection can then be found in the `default_collection` directory. Python CLI has not been tested.

## How it works

...

## Constructing reference sketches

Reference sketches can be rapidly constructed and prepared for use with `Sketchy`. Custom sketches are useful for prediction on species currently not offered in the default collection, lineage sub-sketches of a species, or local genome collections, such as from  healthcare providers or surveillance programs that are not publicly available. All that is needed is a set of high-quality assemblies and their associated genotypes. Ultimately, genome and feature representation in the database should be considered carefully, as they define the genomic neighbors that can be typed with `Sketchy`. 

### Genome assemblies and sketch construction

Assemblies should be of sufficient quality for genotyping and can produced e.g. with classic tools from the [`Torstyverse`](https://github.com/tseemann) like [`Shovill`](https://github.com/tseemann/shovill) or from large-scale public archive surveillance pipelines like [`Pathfinder`](https://github.com/pf-core). 

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

See the [`Nextflow`](#nextflow) section for parallel sketch building and the [`Benchmarks`](#benchmarks) section for guidance on selecting an appropriate sketch and k-mer size for `Sketchy`. 

### Genotype features and index preparation

Genotypes associated with each genome in the reference sketch should be in ta tab-delmited table with a column containing the file name identifiers used in the construction of the reference sketch (for example `uuid`) and it's asosciated genotypes with headers:

```
uuid        st      sccmec  pvl 
DRR083589   st772   v       +
DRR083590   st93    iv      +
DRR119226   st59    v       +
DRR119227   st80    iv      +
DRR128207   st772   -       +
DRR128208   st90    iv      +
```
