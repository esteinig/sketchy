""" ENA assembly collection methods (Blackwell et al. 2021) """

import ijson
import pandas
from pathlib import Path

def parse_metadata(json_file: Path):

    """ Parse species, lineage and quality data of genomes in the ENA collection

    json_file: path to the ENA assembly metadata JSON

    --> retains genomes with full MLST (excluding "-") 

    * filename: Json1_ENA_metadata
    * file address: https://figshare.com/ndownloader/files/26578377
    * data citation: 
        Blackwell, Grace; Hunt, Martin; Malone, Kerri; Lima, Leandro; Horesh, Gal; T. F. Alako, Blaise; et al. (2021): 
        Additional material for "Exploring bacterial diversity via a curated and searchable snapshot of archived DNA sequences". 
        figshare. Dataset. https://doi.org/10.6084/m9.figshare.14061752.v1 
    
    """

    data = []
    with json_file.open() as file:
            parser = ijson.parse(file)
            i = 0
            n = 0
            current_genome = ""
            for prefix, event, value in parser:
                if event == 'start_map' and "." not in prefix:
                    current_genome = prefix
                    bracken = ("", 0.0)
                    completeness = None
                    contamination = None
                    mlst = None
                    mlst_species = None
                    keep = False
                    heterogeneity = None
                if "checkm_results.Completeness" in prefix:
                    completeness = value
                if "checkm_results.Contamination" in prefix:
                    contamination = value
                if "checkm_results.Strain_heterogeneity" in prefix:
                    heterogeneity = value
                if "mlst_results.species" in prefix:
                    mlst_species = value
                if "mlst_results.mlst" in prefix:
                    mlst = value
                    if mlst != "-":
                        keep = True
                if "bracken" in prefix and event == "string" and value not in ("NA", "-"):
                    species = prefix.replace(f"{current_genome}.bracken.", "")
                    abundance = float(value)
                    if abundance > bracken[1]:
                        bracken = (species, abundance)
                if event == 'end_map' and "." not in prefix:
                    if keep:
                        data.append([
                            current_genome, bracken[0], bracken[1], float(completeness), float(contamination), float(heterogeneity), mlst, mlst_species
                        ])
                        n += 1
                    # Reset variables to make sure none are
                    # accidentally used in the next entry
                    bracken = ("", 0.0)
                    completeness = None
                    contamination = None
                    mlst = None
                    mlst_species = None
                    keep = False
                    heterogeneity = None
                # if i == 100000:
                #     continue  # test break
                i += 1

    df = pandas.DataFrame(
        data, columns=["accession", "bracken_species", "bracken_abundance", "completeness", "contamination", "heterogeneity", "mlst", "mlst_species"]
    )
    df.to_csv("meta.tsv", sep="\t", index=False)

def clean_metadata(meta_file: Path, assembly_paths: Path):

    """ Obtain a clean subset of the assembled genomes filtered by assembly quality
    
    * filter genomes: contamination > 1. && completeness < 99. && heterogeneity > 0.1
    * species genome counts, separated into total and > 100
    * FTP paths to species assemblies on the EMBL ENA server

    EMBL EBI address: http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/sampleid_assembly_paths.txt
    
    """

    with meta_file.open() as meta_file:
        df = pandas.read_csv(meta_file, header=0, sep='\t')
        contaminated = df[df['contamination'] > 1.]
        fragmented = df[df['completeness'] < 99.]
        heterogenous = df[df['heterogeneity'] > 0.1]

        all_exclude = pandas.concat((contaminated, fragmented, heterogenous))
        unique_exclude = all_exclude.index.unique()

        df_clean = df[~df.index.isin(unique_exclude)]

        species_counts_100 = []
        species_counts_all = []
        for species, data in df_clean.groupby("bracken_species"):
            n = len(data)
            species_counts_all.append([species, n])
            if n >= 100:
                species_counts_100.append([species, n])
        df_clean.to_csv("meta_clean.tsv", index=False, sep="\t")


        df_species_all = pandas.DataFrame(species_counts_all, columns=["species", "n"]).sort_values("n")
        df_species_100 = pandas.DataFrame(species_counts_100, columns=["species", "n"]).sort_values("n")

        df_species_all.to_csv("species_counts.tsv", index=False, sep="\t")
        df_species_100.to_csv("species_counts_100.tsv", index=False, sep="\t")

        assembly_paths = pandas.read_csv(assembly_paths, sep='\t', header=None, names=['id', 'path'])

        spec_paths = []
        specs = []
        for species, species_df in df_clean.groupby("bracken_species"):
            ftp_paths = assembly_paths[assembly_paths['id'].isin(species_df['accession'])]
            ftp_paths['path'] = [p.replace("/ebi/ftp", "http://ftp.ebi.ac.uk") for p in ftp_paths["path"]]

            species_paths = ftp_paths.drop(columns="id")
            name = f"{species.lower().replace(' ', '_')}"
            species_paths['species'] = [name for _ in species_paths.iterrows()]
            species_df['species'] = [name for _ in species_df.iterrows()]
            spec_paths.append(species_paths)
            specs.append(species_df)

        species = pandas.concat(specs)
        species.to_csv(f"species_data.tsv", sep='\t', index=False)

        species_paths = pandas.concat(spec_paths)
        species_paths.to_csv(f"species_ftp.tsv", sep='\t', index=False, header=False)