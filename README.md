# Sketchy <a href='https://github.com/esteinig'><img src='img/logo.png' align="right" height="210" /></a>

![](https://img.shields.io/badge/version-alpha-red.svg)
![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![](https://img.shields.io/badge/docs-github-green.svg)
![](https://img.shields.io/badge/BioRxiv-v1-orange.svg)

Real-time lineage matching and genotyping from uncorrected nanopore reads

### Overview

**`v0.3-alpha5: public test build, conda install`**

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

---

#### :eyeglasses: `sketchy plot`

Sketchy plot handles the raw output from the prediction and generates a ranked hitmap (by top ranking sum of shared hashes) colored by lineage, and optionally genotype (`-g`) or antimicrobial resistance profiles (`-r`). *K. pneumoniae* does currently not support resistance profiling. Output is a dual plot in the file `{prefix}.{format}` such as `sketchy.png`, where a limit to the reads o nthe `x-axis` can be given by `--limit`:

`sketchy plot -d test_prediction.tsv -p sketchy -f png -g -r --limit 500`

The plot task also generates a plot of the total sum of shared hashes aggregated at each read by lineage, or genotype / resistance profile, which serves as a means of identifying the most frequent value for the trait (ranked in legend), but also serves as comparison tool between the top most common predictions (`--top`) across the window outlined by the hitmap. Colors (`--color`) can be `brewer` palette names such as `PuGn` or a comma delimited list of `brewer` palette names if for example genotype (`-g`) is activated:


`sketchy plot -d test_prediction.tsv -p sketchy -f png -g --top 5 --color Blues_r,Greens_r`

When the breakpoint `-b` option is activated the task attempts to determine a breakpoint on the most frequent trait where the sum of sum of shared hashes (2nd plot) is stable for `--stable` amount of reads. This threshold by default is set to 500, but may need to be adjusted for species liek *M. tuberculosis*. This option will also output a file `{prefix}.bp.tsv` which writes the breakpoints for later parsing, for example i nthe bootstrap Nextflow


`sketchy plot -d test_prediction.tsv -p sketchy -f png -b --stable 500`

---

#### :umbrella: `nextflow sketchy.nf`

Bootstrapping workflow in Sketchy.
