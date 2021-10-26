![sketchy logo](images/logo.png)


`Sketchy` is a lineage calling and genotyping tool based on the heuristic principle of genomic neighbor typing developed by [Karel Břinda and colleagues (2020)](https://www.biorxiv.org/content/10.1101/403204v2). It queries species-wide ('hypothesis-agnostic') reference sketches using MinHash and infers associated genotypes based on the closest match, including multi-locus sequence types, susceptibility profiles, virulence factors or other genome-associated features provided by the user. Unlike the original implementation, `Sketchy` does not use phylogenetic trees (which has some downsides) and is easier to modify for users, for example to construct reference sketches for local genome collections, or species not implemented in the current release.

## Installation

<div class="termy">

```console
$ cargo install sketchy

---> 100%
```

</div>

## Help

Show the subcommands available for Skethy

<div class="termy">

```console
$ sketchy

sketchy 0.6.0
Predict genotypes based on cumulative shared hashes

USAGE:
    sketchy <SUBCOMMAND>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    check      
    help       Prints this message or the help of the given subcommand(s)
    info       
    predict    
    shared     
    sketch
```

</div>

## Example

For demonstration we will build a high-resolution reference sketch of *Pseudomonas aeruginosa* isolates from the ENA reference assembly collection. For a full guide on how to simulate test data sets with `Sketchy` see the [`Walkthrough`]('walkthrough/building_reference_sketches.md) section.
