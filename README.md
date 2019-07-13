# Sketchy <a href='https://github.com/esteinig'><img src='img/logo.png' align="right" height="250" /></a>

Real-time lineage matching and genotyping from uncorrected nanopore reads

![](https://img.shields.io/badge/version-alpha-red.svg)
![](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
![](https://img.shields.io/badge/docs-latest-green.svg)
![](https://img.shields.io/badge/BioRxiv-prep-green.svg)

## Overview

**`v0.2-alpha: internal pre-release, no tests`**

`Sketchy` is an online lineage matching algorithm for real-time genotyping and susceptibility prediction in bacterial pathogens using nanopore sequencing platforms. `Sketchy` uses MinHash distances to match reads to a representative population genomic sketch of the target organism. Because it relies on well characterized, high-quality genomes from the public domain, it is currently restricted to common pathogens for which we have sufficient data that can be processed with the [`pf-core/pf-survey`](https://github.com/pf-core) Nextflow pipeline, which at this stage supports *Staphylococcus aureus*, *Klebsiella pneumoniae* and *Mycobacterium tuberculosis* specific typing. *Burkholderia pseudomallei* with antibiotic resistance genotype using `ArDAP` and other candidates are planned for beta release.

Preprint coming soon.

## Install

```
git clone https://github.com/esteinig/sketchy
pip install sketchy
```

## Command line interface

##### :eyeglasses: `sketchy sketch`

To build a MinHASH sketch with MASH from a collection of `FASTA` assemblies with associated results from genotyping or pohenotype data:

```
Usage: sketchy sketch [OPTIONS]

  Create a MinHash sketch for matching with MASH

Options:
  -d, --data TEXT     Input data file to assign lineages and genotypes
  -k, --kmer INTEGER  K-mer length to create sketch with in MASH
  -s, --size INTEGER  Sketch size in MASH
  -o, --outdir TEXT   Output directory for sketch files
  -p, --prefix TEXT   Prefix for data sketch with Sketchy
  -c, --copy          Copy assembly files
  --delimiter TEXT    Delimiter for data file
  --help              Show this message and exit
```

---

##### :briefcase: `sketchy predict`

Apply the sketchy algorithm to a cumulatively sliced `FASTQ` files with basecalled nanopore reads. Warning: this may generate a collection of temporary files with cumulative reads at the moment, depending on the `--reads` number of reads to type. Do not type the whole run at the moment, at most predictions are definitely stable around 200 - 300 reads.

```
Usage: sketchy predict [OPTIONS]

  Online lineage matching from uncorrected nanopore reads

Options:
  -f, --fastq TEXT     Input Fastq file to predict from
  -s, --sketch TEXT    MASH sketch to query
  -d, --data TEXT      Index data file to pull genotypes and other data from
  -t, --tmp TEXT       Temporary dir for slicing Fastq file
  -c, --cores INTEGER  Number of processors for `mash dist`
  -h, --header         Print header to online mode STDOUT.
  -r, --reads INTEGER  Number of reads to type.
  --single             Single read analysis, does not compute score
  --dist               Use best hits from min hash distance
  -o, --output TEXT    Output CSV file for --single
  --help               Show this message and exit.
```

---

##### :closed_umbrella: `sketchy sort`

Apply to sort a `FASTQ` file by the start time tag in the header in files called by Albacore.

```
Usage: sketchy sort [OPTIONS]

  Sort basecalled reads by start time (Albacore)

Options:
  -f, --fastq TEXT   Input FASTQ file to slice and predict
  -o, --output TEXT  Sorted by start datetime FASTQ output file
  -s, --shuffle      Shuffle start datetimes in FASTQ
  --help             Show this message and exit.
```
