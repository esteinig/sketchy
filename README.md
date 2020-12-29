# sketchy <a href='https://github.com/esteinig'><img src='docs/logo.png' align="right" height="250" /></a>

![](https://img.shields.io/badge/lang-rust-black.svg)
![](https://img.shields.io/badge/version-0.5.0-purple.svg)
![](https://img.shields.io/badge/biorxiv-1.0-blue.svg)

Real-time lineage hashing and genotyping of bacterial pathogens from uncorrected nanopore reads using genomic neighbor typing with `Mash`

## Overview

**`v0.5.0: preprint`**

`Sketchy` is a lineage calling and genotyping platform based on the heuristic principle of genomic neighbor typing developed by [Karel Břinda and colleagues (2020)](https://www.biorxiv.org/content/10.1101/403204v2). `Sketchy` wraps `mash screen` for completed sequence runs and a streaming version of `mash dist` (the sum of shared hashes) for real-time analysis. It queries species-wide, lineage-resolved reference sketches of bacterial whole genome assemblies and infers their associated genotypes based on the closest reference matches, including multi-locus sequence types, susceptibility profiles, virulence factors or species-specific markers. Precomputed genotype features and automated updates of species databases can be found in the corresponding pathogen reference sections. 

Species we have validated using ONT sequence reads with matching Illumina data:

* [*Staphylococcus aureus*](docs/saureus.md) (n = 142)
* [*Klebsiella pneumoniae*](docs/kpneumo.md) (n = 120)

**2021**: sketch updates will occur every month and include new quality controlled genomes from the European Nucleotide Archive

Please see our preprint for guidance on the limitations of `Sketchy`.

- [Install](#install)
- [Usage](#sketchy-usage)
- [Reference sketches](#reference-sketches)
- [Constructing reference sketches](#reference-sketches)
- [Tasks and parameters](#tasks)
- [Citing](#citing)

## Install

Sketchy implements a Rust command-line interface (`sketchy-rs`) for computation and evaluation on read streams and a Python command-line interface (`sketchy`) for evaluation plots and database utilities.

**`Cargo`**

Rust client only where the `compute` subtask requires `Mash` in host `$PATH`

```bash
cargo install sketchy-rs
```

On Linux systems one should be able to install `Mash` conveniently e.g.

```bash
sudo apt-get install mash
```

**`Singularity`**

Singularity container is preloaded with default reference sketches at the root path `/sketchy`

```sh
singularity pull docker://esteinig/sketchy:latest
```

**`Docker`**

Docker image is available, also comes preloaded with the default reference sketches

```sh
docker pull esteinig/sketchy:latest
```

**`BioConda`**

`Sketchy` is available on `BioConda` (thanks to [@mbhall88](https://github.com/mbhall88))

```
conda install -c bioconda -c conda-forge sketchy=0.5.0
```

You can also manually install the latest commits into an environment like this:

```sh
conda install -c bioconda -c conda-forge \
  mash=2.2.2 psutil pysam rust nanoq
cargo install sketchy-rs
git clone https://github.com/esteinig/sketchy
pip install ./sketchy
sketchy pull
sketchy list
```

## Reference sketches

If no container is used, pull default species sketches into local storage before first use:

```
sketchy pull
```

You can set the path to which the sketches are downloaded e.g. to the default `--path ~/.sketchy`. 

```
sketchy pull --path ~/.sketchy
```
Local sketches and template names can be viewed with:

```
sketchy list
```

Set the environment variable `$SKETCHY_PATH` to a custom sketch directory for tasks to discover databases automatically.

## Sketchy usage

See the `Tasks and Parameters` section for details on all tasks and settings available in `Sketchy`. Reads are expected to belong to the species of the selected reference sketch. For an evaluation of genomic neighbor typing using `mash screen` and a comparison to the online implemention of `mash dist` on comprehensive and dereplicated strain-level sketches, please see the preprint. 

Setup environment variable to databases (or specify database paths in `--db` arguments)

```bash
SKETCHY_PATH=$HOME/.sketchy
```

### Screening function

`Sketchy` primarily uses a screening of the reference sketch containment in the provided read set wrapping `Mash`. I tend to use this function for quick and easy genomic neighbor type screening on many isolates. Screening with `Mash` uses the winner-takes-all strategy and `Sketchy` then simply links the best match with the genotype data provided with the reference sketches. 

In the Python client:

```
sketchy screen --fastx test.fastq --db saureus --limit 10 --pretty
```

In the Rust client:

```
sketchy-rs screen -f test.fastq -d saureus -l 10 -p
```

Please cite the following when using `sketchy screen`:

* Ondov et al. (2016) - `Mash`
* Ondov et al. (2019) - `Mash Screen`

### Streaming function

Streaming genomic neighbor typing heuristic that implements `mash dist` and computes the sum of shared hashes against the reference sketch. Because streaming is slower than screening for completed sequence runs, I tend to use this more in cases where few reads are available or when streaming is actually required (not that often). In some edge cases the streaming utility can be quite useful - for instance, we confirmed a *S. aureus* re-infection of the same strain in a cystic fbrosis patient from < 27 reads and diagnostic plots, which was not possible with the `screen` implementation.

`Sketchy's` streaming algorithm can be run through a wrapper in the Python client which is only suitable for completed read files (`sketchy run`). Read streams and online sequencing runs should be served with the Rust client (see below). 

```bash
sketchy run --fastq test.fq --limit 1000 --ranks 10 --outdir test
```

Streaming is primarily bottlenecked by sketch queries of each read against the reference sketch, which means that prediction speeds are usually fast on smaller sketches (e.g. 10,000 genomes, ~ 100 reads/second) but for large sketches and analyses over tens of thousands of reads total runtime can be excruciating. Fortunately, we generally do not need that many reads to make predictions. When using species-wise reference sketches with tens of thousands of genomes on large read sets use `head` or `--limit` options in the command line clients to predict on the first few thousands reads (`sketchy run --limit 1000` or `cat test.fq | head -4000 | sketchy-rs compute`) which should be sufficient for initial analysis. Smaller reference sketches by lineage or created from local collections should be sufficiently fast for online prediction on MinION / Flongle / GridION.


The `Rust` command line interface implements two subtasks: `sketchy-rs stream` which computes sum of shared hashes and ranked sums of shared hashes by genotypes, and `predict` which uses the output to predict the genotype profile. 

 ```bash
 cat test.fq | sketchy-rs stream -d saureus -t 4 > sssh.tsv
 ```
 
`Predict` then uses the ranked scores at each read to infer the genotype profiles:

```bash
cat test.ssh.tsv | sketchy-rs predict -d saureus > predict.tsv
```

`Stream` and `predict` read from `/dev/stdin` so they can be piped:

```bash
cat test.fq | sketchy-rs stream -d saureus | sketchy-rs predict -d saureus > predict.tsv
```

Diagnostic plots and evaluation summaries are handled in the Python client and accessed via the `sketchy plot` task:

```
sketchy plot \
    --sssh sssh.tsv \
    --db saureus \
    --palette YnGnBu \
    --prefix test \
    --format png
```

Please cite the following when using `sketchy stream`:

* Ondov et al. (2016) - `Mash`

### Online streaming analysis

In a live sequencing run, `Sketchy` can be set to observe a directory (e.g. `fastq_pass` from live basecalling) in order to stream reads into the Rust client. A watcher waits for the `fastq` file to be completed before piping the filename to `/dev/stdout` and into the streaming client:

```
sketchy online watch -d /live/fastq | cat - | sketchy-rs stream --db saureus > test.sssh.tsv
```

### Android mobile phones

To set up the Rust client on Android mobile phones, the following can be done quite easily:

1. Install the [`UserLAnd`](https://github.com/CypherpunkArmory/UserLAnd) app 
2. Setup an `Ubuntu` image
3. Run the following in the terminal

```bash
sudo apt-get update && sudo apt-get install curl mash build-essential

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
cargo install sketchy-rs

wget https://storage.googleapis.com/sketchy-sketch/saureus.min.tar.gz \
  && tar -xvzf saureus.min.tar.gz
wget https://storage.googleapis.com/sketchy-sketch/mobile_test.fq

cat mobile_test.fq | sketchy-rs stream --db ./saureus.min
```

### Nextflow pipeline

TBD.
 

## Reference sketches

Reference sketches can be constructed and prepared for use with `Sketchy`. Custom sketches are useful for prediction on species currently not offered in the default collection, sub-sketches of a species, or local genome collections, such as from healthcare providers or surveillance programs not publicly accessible. All that is required is a set of high-quality assemblies and their associated genotypes. Ultimately, genome and feature representation in the database should be considered carefully, as they define the genomic neighbors that can be typed with `Sketchy`. 

Custom sketches must include a:

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

### Genome assemblies and sketch construction

Given a set of high-quality assemblies in the format `{id}.fasta`:


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
mash sketch -s 1000 -k 15 -o ref *.fasta
```

See the [`Nextflow`](#nextflow) section for parallel sketch building and the [`Benchmarks`](#benchmarks) section for guidance on selecting an appropriate sketch and k-mer size for `Sketchy`. 

### Genotype features and index preparation

Genotypes associated with each genome in the reference sketch should be in a tab-delimited table (`genotypes.tsv`) with appropriate headers. We will remove the `id` column in the next step when we translate columns into genotyping indices for `Sketchy`.

```
id          st      sccmec  pvl 
DRR083589   st772   v       +
DRR083590   st93    iv      +
DRR119226   st59    v       +
DRR119227   st80    iv      +
DRR128207   st772   -       +
DRR128208   st90    iv      +
```

To generate the `Sketchy` reference genotypes in `sketchy genotypes create`:

```
sketchy genotypes create -i genotypes.tsv -s ref.msh --prefix ref --drop uuid
```

This will create a directory with the following files

```
ref.msh   # sketch
ref.tsv   # genotypes
ref.idx   # index
ref.json  # key 
```

## How the streaming algorithm works

`Sketchy` stream computes across three stages and two simple scores: the first is the sum of shared hashes, where it keeps a cumulative sum of shared hashes (`ssh`) computed against each index in the reference sketch for each consecutive read. `Mash` queries take the majority of the compute while `Sketchy` siphons off the output: at each read, the indexed scores are ranked and the highest ranking scores are recorded to reduce excessive read-wise output of the reference sketch queries and retain salient hits (default: `--ranks 10`). Raw sum of shared hashes scores can be output for debugging using the `--ssh` scores flag.

In the second stage, an evaluation score for genotypes features is computed by aggregating the sum of ranked sums of shared hashes (`sssh`) for each genotype feature in the associated index independently (e.g. *mecA* presence or SCC*mec* type). Evaluations are plotted for visual confirmation, along with a preference score adopted from [Brinda and colleagues](https://www.biorxiv.org/content/10.1101/403204v2) that indicates the degree of confidence in the best prediction over the second-best prediction.

In the third stage, the sum of ranked sums of shared hashes scores are ranked and translated to genotypes from the database. As each genotype feature is evaluated independently, ranked genotype predictions are combined at each rank to provide a full genotype. For example, it may be that two ranks of sequence types are predicted (e.g. ST1 and ST8), but only one rank of penicillin resistance (R). In this case, the first rank genotype is unambigious (e.g. ST1-R), but penicillin resistance is adopted in the second rank genotype (e.g. ST8-R)

### Sketchy evaluation outputs

Sketchy produces a directory `--output` with the intermediary pipeline data files (`prefix.ssh.tsv` and `prefix.sssh.tsv`). For evaluation and prediction output, the primary data file is `prefix.data.tsv` which shows the final prediction for each genomic feature, the determined stability breakpoints in reads (`0` or `-1` in < v0.4.4 means that a breakpoint could not be called either because predictions were not stable or the chosen stable breakpoint was smaller than the evaluated reads) and the median preference score over the evaluated reads:

```
feature         prediction      stability       preference
mlst            ST93            23              0.66666667
meca            MRSA            17              0.39714868
pvl             PVL+            23              1.0
scc             SCCmec-IV       23              0.39745917
clindamycin     S               1               0.80033841
rifampicin      S               1               1.0
ciprofloxacin   S               1               1.0
vancomycin      S               1               1.0
tetracycline    S               2               1.0
```

Debugging plots for evaluation are the more salient outputs - the following images show a different sample prediction than the outputs above. Each row in the image corresponds to one genomic feature prediction, which is listed in the middle legend together with the top five alternative feature value predictions. Each feature value prediction corresponds to a color, where darker colors represent the highest-ranking and therefore most likely predictions.

<a href='https://github.com/esteinig'><img src='docs/example_saureus_1.png' align="center" height="500" /></a>

**What's going on here?**

In the heatmap, the highest-ranking (descending) raw sum of shared hashes queries against the database sketch are shown and colored. Gray colors represent feature values not in the ultimate highest-ranking five and demonstrate uncertainty in the initial predictions. On the other hand, homogenous color represents certainty in the prediction which may increase as the scores are updated.

In the middle plot, the ranked sum of shared hashes (`ssh`) are evaluated by aggregating the sum of their ranked sum of shared hashes (`sssh`) by feature value, from which stability breakpoints are calculated (vertical lines) i.e. where the highest scoring feature value remains the highest scoring for `--stable` reads. In the example this defaults to 1000 reads, so no breakpoints were detected (set to `0`) as the prediction was limited to 1000 reads total; breakpoints are included in the `prefix.data.tsv` output file. Legend items and colors are ordered according to rank; a straight, uncontested line for a dominant feature value score indicates certainty the same as homogenous color in the heatmap.

In the plot on the right, the preference score from [Brinda and colleagues](https://www.biorxiv.org/content/10.1101/403204v2) is computed on the sum of ranked sums of shared hashes (`sssh`) scores from the middle plot. As in the original a threshold of `p = 0.6` (horizontal line) indicates when a prediction should be trusted and when it should not. Note that the preference is always computed on the feature value with the highest score over the feature value with the second highest score, regardless of whether it is the right prediction. In fact, the score is susceptible to 'switches' in predictions, especially using lower resolution sketches, where a prediction is updated and flips to another more likely prediction as more evidence is gathered. 

<a href='https://github.com/esteinig'><img src='docs/example_saureus_2.png' align="center" height="500" /></a>

In this example, the same data from the Bengal Bay clone is run on the lower resolution reference sketch `saureus_15_1000` instead of `saureus_15_10000`. Incorrect sequence type ST12 is called for about 300 reads before making a switch to the correct sequence type ST772. This is reflected in the heatmap by distinct color blocks, but lower-resolution also trades-off prediction speed with larger more accurate sketches. In the higher resolution sketch above, the sequence type is called almost immediately and initial uncertainty is lower, as indicated by less gray coloring in the heatmap on the initial reference sketch queries.
