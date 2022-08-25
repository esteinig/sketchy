import ijson
import pandas

DATA_MLST_COMPLETE = False

if DATA_MLST_COMPLETE:
    data = []
    with open("assembly.ena.json") as file:
            parser = ijson.parse(file)
            i = 0
            n = 0
            current_genome = ""
            for prefix, event, value in parser:
                if event == 'start_map' and "." not in prefix:
                    current_genome = prefix
                    data_current_genome = []
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
                        print('KEEP', current_genome, n)
                    # Reset variables to make sure none are
                    # accidentally used in the next entry
                    bracken = ("", 0.0)
                    completeness = None
                    contamination = None
                    mlst = None
                    mlst_species = None
                    keep = False
                    heterogeneity = None
                if i == 100000:
                    continue  # test break
                i += 1

    df = pandas.DataFrame(
        data, columns=["accession", "bracken_species", "bracken_abundance", "completeness", "contamination", "heterogeneity", "mlst", "mlst_species"]
    )
    df.to_csv("meta.tsv", sep="\t", index=False)

PROCESS_SPECIES_SKETCHES=False

SPEC = ("Staphylococcus aureus", "Streptococcus pneumoniae", "Neisseria meningitidis", "Klebsiella pneumoniae", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa")

if PROCESS_SPECIES_SKETCHES:
    with open("meta.tsv", "r") as meta_file:
        df = pandas.read_csv(meta_file, header=0, sep='\t')
        contaminated = df[df['contamination'] > 1.]
        fragmented = df[df['completeness'] < 99.]
        heterogenous = df[df['heterogeneity'] > 0.1]

        all_exclude = pandas.concat((contaminated, fragmented, heterogenous))
        unique_exclude = all_exclude.index.unique()

        df_clean = df[~df.index.isin(unique_exclude)]

        print(df)
        print(df_clean)

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
        
        print(df_species_all.n.sum())
        print(df_species_100.n.sum())

        assembly_paths = pandas.read_csv("assembly_paths.txt", sep='\t', header=None, names=['id', 'path'])


        for species in SPEC:
            species_df = df_clean[df_clean["bracken_species"] == species]
            ftp_paths = assembly_paths[assembly_paths['id'].isin(species_df['accession'])]
            ftp_paths['path'] = [p.replace("/ebi/ftp", "http://ftp.ebi.ac.uk") for p in ftp_paths["path"]]

            species_paths = ftp_paths.drop(columns="id")
            name = f"{species.lower().replace(' ', '_')}"
            species_df.to_csv(f"{name}.tsv", sep='\t', index=False)
            species_paths.to_csv(f"{name}.ftp.tsv", sep='\t', index=False, header=False)

SUBSET_SAUREUS=True

if SUBSET_SAUREUS:

    accessions = pandas.read_csv("run_query.tsv", header=0, sep='\t')
    saureus = pandas.read_csv("saureus.tsv", header=0, sep='\t')
    saureus_ena = pandas.read_csv("staphylococcus_aureus.tsv", header=0, sep='\t')
    print(saureus_ena)

    common = accessions[accessions['sample_accession'].isin(saureus_ena['accession'])]

    clean_genotypes = saureus[saureus['id'].isin(common['run_accession'])]
    clean_genotypes_merged = clean_genotypes.merge(common, left_on="id", right_on="run_accession")

    sketch_genotypes = clean_genotypes_merged.iloc[:, :16]
    sketch_genotypes.insert(loc=0, column='id', value=clean_genotypes_merged['sample_accession'].tolist())

    clean_genotypes_merged.to_csv("saureus.ena2.tsv", sep='\t', index=False)
    clean_accessions = clean_genotypes_merged[['id']]
    clean_accessions.to_csv("saureus.ena2.ids.tsv", sep='\t', index=False)

    sketch_genotypes = sketch_genotypes[~sketch_genotypes.duplicated()]
    sketch_genotypes = sketch_genotypes[~sketch_genotypes['id'].duplicated()]
    sketch_genotypes.to_csv("saureus.ena2.sketch.tsv", sep='\t', index=False)


ID_SAUREUS=False
SIMSPEC = ("staphylococcus_aureus", "streptococcus_pneumoniae", "neisseria_meningitidis", "klebsiella_pneumoniae", "neisseria_gonorrhoeae", "pseudomonas_aeruginosa")
if ID_SAUREUS:

    for species in SIMSPEC:
        simref = pandas.read_csv(f"species/simulated/{species}.sim.tsv", header=None, sep='\t', names=["run"])
        simref_names = simref["run"].tolist()

        species_mlst = pandas.read_csv(f"species/{species}.tsv", header=0, sep='\t')

        species_simref = species_mlst[species_mlst["accession"].isin(simref_names)]

        print(len(simref_names), len(species_simref))

        species_simref[["accession", "mlst"]].to_csv(f"species/simulated/{species}.simref.tsv", index=False, sep="\t", header=True)

