# Extracting validation data from the Blackwell collection

This process extracts the MLST validation data from the Blackwell collection as outlined in the manuscript.

1. Download the meta-data JSON (https://figshare.com/ndownloader/files/26578377) 
2. Run the `species.py` application (`python species.py --help`) on the JSON file, this will produce a summary of the meta-data including MLST (`species_data.tsv`), a file for all FTP addresses (`species_ftp.tsv`) and count data for all species (Bracken, `species_counts.tsv`)
3. Download the assemblies from the FTP addresses for a species of interest, these can be used for sketch construction. Genotype files can be constructed from subsets of the `species_data.tsv` file for the assemblies included in the sketch. 
4. You may want to run a check on the order and congruence between sketch and genotype file. Commands for sketch construction are outlined in the `local sketches` and `genotype files` sections of the [documentation](https://github.com/esteinig/sketchy/blob/master/docs/index.md).
