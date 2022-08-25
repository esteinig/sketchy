## Overview

`Sketchy` is a nanopore lineage calling and genotyping tool based on the heuristic principle of [genomic neighbor typing (Břinda et al. 2020)](https://www.biorxiv.org/content/10.1101/403204v2). `Sketchy`  queries species-wide (hypothesis-agnostic) reference sketches using MinHash methods to infer genotypes based on the closest reference match. Reference databases and genotypes, such has multi-locus sequence types, susceptibility profiles, or virulence factors, are configurable by users.

## Install

Rust client.

```
$ cargo install sketchy
```

Rust client and Python utility client.


```
$ conda install -c bioconda sketchy
```

[Release binaries](releases.md) for Linux and MacOS are available.


## Subcommands

Show subcommands available in the Rust client.

```
$ sketchy --help
```

Show subcommands available in the Python client.

```
$ sketchy-utils --help
```

## Predictions

`Sketchy` predicts the genotype of the genome with the highest shared hashes in the reference sketch. Two modes are available:

1. **Offline (read sets)**: k-mers are hashed and shared hashes computed from a set of reads. Output are the `--top` genome matches and their genotypes. 
2. **Online  (streaming)**: k-mers are hashed for each incoming read, the sum of shared hashes for each genome in the reference sketch is updated and the `--top` genome matches and their genotypes based on the sorted sum of shared hashes at the current read are printed.

Sum of shared hashes based on streaming are slightly less informative than shared hashes based on read sets. Streaming requires less memory, but can be slow for large read sets, especially when deploying large reference sketches.

Predictions require the input sequences (`-i`), a reference sketch (`-r`) and a matching genotype file (`-g`) as described below ([genotype files](#genotype-files)).

### Read sets

Output the best 5 matches against the reference sketch with header.

```console
$ sketchy predict -i seq.fq -r saureus.msh -g saureus.tsv -t 5 -H
```

### Streaming

Output the updated best match against the reference sketch from a stream of reads.

```console
$ cat seq.fq | sketchy predict -r saureus.msh -g saureus.tsv -s
```

## Sketches

### Species sketches

Species sketches are available in the [data repository]() and can be built from local reference assembly collections.

### Local sketches

Reference sketches can be built from any collection of assembled genomes for which associated genotype or phenotype data are available. 

Build a default-resolution (`s = 1000`) reference database from any collection of `fasta`.

```
$ sketchy sketch -i *.fa -k 16 -s 1000 -o ref.msh
```

You can pipe assemblies using `find` into the sketch construction, as wildcard expansions are limited to ~30,000 files.

```
$ find assemblies/ -name "*.fa" | sketchy sketch -k 16 -s 1000 -o ref.msh
```

List sketch parameters.

```
$ sketchy info -i ref.msh -p 
```

### Genotype files

Prediction requires a **tab-delimited** genotype index **in the same order and of the same length** as the reference sketch. Names in the genotype index (first column) are the file name of the input files in the sketch

```
name    mlst    tetracyline penicillin methicilin
ERR129347.fa    ST82    R   R   S
ERR121347.fa    ST93    S   S   S
```

List the order of genomes in the sketch, their length (bp) and an estimate of cardinality (bp)

```
$ sketchy info -i ref.msh
```

Check if the genotype file contains the correct number and order of genomes as the sketch.

```
$ sketchy check -r ref.msh -g ref.tsv
```

This will outpout `ok` to `stdout` if the check is completed successfully and fail with an error message otherwise.

## Sketch validation

We conducted simulations using `badread` of the following species:

* *Neisseria spp.*
* *Streptococcus pneumoniae*
* *Klebsiella pneumoniae*
* *Staphylococcus aureus*
* *Pseudomonas aeruginosa*

Cross-validation sketches (using subsampling of reference assemblies) and read data are available in the [data repository]().

## Other


### Shared hashes

Given two assembled genome sequences, create a sketch at default k-mer size of `k = 16` and sketch size of `s = 1000`. `Mash` configuration can be replicated by setting `--seed 42`.

Sketch two genome assemblies with identical settings.

```
$ sketchy sketch -i genome1.fa -o genome1.msh
$ sketchy sketch -i genome2.fa -o genome2.msh 
```

Compute shared hashes between the reference and query genomes.

```
$ sketchy shared -r genome1.msh -q genome2.msh

> genome1.fa genome2.fa 360
```


If multiple sketches are available compute pairwise shared hashes.

```
$ sketchy sketch -i genome1.fa genome2.fa -o multi.msh
$ sketchy shared -r multi.msh -q genome2.msh

> genome1.fa genome2.fa 360
> genome2.fa genome2.fa 1000
```
