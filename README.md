# Sketchy <a href='https://github.com/esteinig'><img src='img/logo.png' align="right" height="210" /></a>

![](https://img.shields.io/badge/version-alpha-red.svg)
![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![](https://img.shields.io/badge/docs-github-green.svg)
![](https://img.shields.io/badge/BioRxiv-v1-orange.svg)

Real-time lineage matching and genotyping from uncorrected nanopore reads

### Overview

**`v0.3-alpha3: public test build`**

`Sketchy` is an online lineage matching algorithm for real-time genotyping and susceptibility prediction in bacterial pathogens using nanopore sequencing platforms. Currently supported species are *Staphylococcus aureus*,  *Klebsiella pneumoniae* and *Mycobacterium tuberculosis*.

### Install
---

* :snake: `conda install -c bioconda -c conda-forge -c esteinig sketchy`
* :whale: `docker pull esteinig/sketchy`

### Usage
---

#### :briefcase: `sketchy predict`

Main interface for prediction on uncorrected nanopore reads. Sketches (`-s`) available are: *S. aureus* (`mrsa`), *K. pneumoniae* (`kleb`) and *M. tuberculosis* (`tb`)

`sketchy predict --help`

Completed sequence read file (`.fastq`) - predict on first 1000 reads (default) and compute the sum of shared hashes post-hoc, 8 processors, using the *S. aureus* sketch:

`sketchy predict -f test.fq.gz -t 8 -s mrsa -o test_prediction`

This produces the data file `test_prediction.tsv` which is the input for `sketchy plot`

#### :eyeglasses: `sketchy plot`


This produces the 
