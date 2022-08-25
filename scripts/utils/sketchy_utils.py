
from math import nan
import typer
import random
import pandas
import shutil

import seaborn as sns
from matplotlib import pyplot as plt
from statistics import mean, stdev
from typing import Optional
from pathlib import Path

app = typer.Typer(add_completion=False)


SPECIES_DICT = {
    "sa": "S. aureus",
    "nm": "N. meningitidis",
    "ng": "N. gonorrhoeae",
    "pa": "P. aeruginosa",
    "kp": "K. pneumoniae",
    "sp": "S. pneumoniae"
}

@app.command()
def sample(
    directory: Path,
    samples: Optional[int] = typer.Option(None, "-n"),
    proportion: Optional[float] = typer.Option(0.01, "-p"),
    extension: Optional[str] = typer.Option(None, "-e"),
    move_to: Optional[Path] = typer.Option(None, "-m"),
    copy_to: Optional[Path] = typer.Option(None, "-c"),
):
    """
    Take a random sample of assembly files for simulation
    """

    files = [f for f in directory.glob(f"*{extension if extension is not None else ''}") if f.is_file()]
    
    if samples:
        n = samples
    elif proportion:
        n = int(proportion*len(files))
    else:
        print("No sampling proportion or number of samples given")
        exit(1)

    sampled = [f for f in random.sample(files, n)]
    
    print(f"Sampled {len(sampled)} / {len(files)}")
    if move_to:
        if not move_to.exists():
            move_to.mkdir(parents=True)
        for f in sampled:
            f.rename(move_to / f.name)
    elif copy_to:
        if not copy_to.exists():
            copy_to.mkdir(parents=True)
        for f in sampled:
            shutil.copy(f, copy_to / f.name)
    else:
        for f in sampled:
            print(f)


@app.command()
def simeval(
    eval_dir: Path,
    truth_dir: Path
):

    species_df = {
        file.stem: pandas.read_csv(file, sep="\t") 
        for file in truth_dir.glob("*.tsv")
    }  # contains the simulated random sample data

    data = []
    for file in eval_dir.glob("*.tsv"):
        identifier = file.stem.split("_")
        species, reads = identifier[0] + "_" + identifier[1], identifier[2]

        print(species, species_df.keys())
        assert species in species_df

        prediction = pandas.read_csv(file, sep=" ", header=None, names=['name', 'mlst']).sort_values('name')
        
        species_data = species_df[species]
        species_truth = species_data[species_data['accession'].isin(prediction.name)].sort_values('accession')

        truth_mlst = species_truth['mlst'].tolist()
        truth_name = species_truth['accession'].tolist()
        predict_mlst = prediction['mlst'].tolist()
        predict_name = prediction['name'].tolist()

        assert sum([1 if name_truth == name_predict else 0 for name_truth, name_predict in zip(truth_name, predict_name)]) == len(species_truth)
        truth_vec = [1 if mlst_truth == mlst_predict else 0 for mlst_truth, mlst_predict in zip(truth_mlst, predict_mlst)]        
        truth_proportion = sum(truth_vec) / len(truth_vec)
        data.append([species, int(reads), truth_proportion])
    
    df = pandas.DataFrame(data, columns=['species', 'reads', 'correct']).sort_values(["species", "reads"])

    
    fig, axes = plt.subplots(
        nrows=1, ncols=1, figsize=(14, 10)
    )

    sns.lineplot(
        data=df, x="reads", y="correct", hue="species", ax=axes, palette="colorblind", marker="o"
    )

    plt.tight_layout()
    fig.savefig(f"eval.png")

@app.command()
def simeval_species(
    eval_dir: Path,
    reads: int = typer.Option(1000, "-r"),
    ext: str = typer.Option("png", "-e"),
    sketch_size: int = typer.Option(1000, "-s"),
):
    
    _, summary_1000 = _get_sketchy_summary(eval_dir)

    pal = sns.color_palette('colorblind') 
    fig, axes = plt.subplots(
        nrows=1, ncols=1, figsize=(14, 10)
    )
    summary_1000 = summary_1000.replace(['simdb', 'refdb'], ['excluded from database', 'included in database'])
    p1 = sns.violinplot(
        data=summary_1000, x="species", y="correct", hue="db", ax=axes, palette=[pal[2], pal[-1]], inner="quartile", split=True, order=[
            "N. gonorrhoeae", "N. meningitidis", "S. pneumoniae", "K. pneumoniae", "S. aureus", "P. aeruginosa"
        ]
    )
    p2 = sns.swarmplot(
        data=summary_1000, x="species", y="correct", hue="db", ax=axes, color="black", dodge=True, order=[
            "N. gonorrhoeae", "N. meningitidis", "S. pneumoniae", "K. pneumoniae", "S. aureus", "P. aeruginosa"
        ]
    )
    handles, labels = axes.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], loc=4, title='Genomes', fontsize='x-large', title_fontsize='16')
    plt.setp(axes.collections, alpha=.7)
    plt.title(f"MLST prediction at {reads} reads [s = {sketch_size}]", size=16, fontweight="bold")
    plt.xlabel("")
    plt.ylim((-0.2, 1.2))
    
    _, xticks = plt.xticks()
    p2.set_xticklabels(xticks, size = 14)
    p1.set_yticklabels([round(f, 2) for f in p1.get_yticks()], size = 14)

    plt.ylabel("Proportion correctly predicted [n = 20, 10 replicates]", size=16)
    plt.tight_layout()
    
    fig.savefig(f"species_sketchy_{reads}.{ext}")

@app.command()
def simeval_method(
    eval_dir: Path
):

    flye_summary = _get_mlst_summary(eval_dir, method='flye')
    medaka_summary = _get_mlst_summary(eval_dir, method='medaka')

    _, summary_1000 = _get_sketchy_summary(eval_dir)
    pal = sns.color_palette('colorblind') 

    krocus_summary = _get_krocus_summary(eval_dir)

    # Plot of mean truth proportion

    fig2, axes2 = plt.subplots(
        nrows=1, ncols=1, figsize=(14, 10)
    )

    summary_1000_sim = summary_1000[summary_1000['db'] == 'simdb']
    summary_1000_sim['method'] = ['sketchy' for _ in summary_1000_sim.iterrows()]

    print(summary_1000_sim)
    methods = pandas.concat([summary_1000_sim, krocus_summary, flye_summary, medaka_summary])
    print(methods)

    p = sns.barplot(
        data=methods, x="method", y="correct", hue="species", ci="sd", ax=axes2, palette=[pal[0], pal[1], pal[2], pal[3], pal[8], pal[9]], alpha=0.92, hue_order=[
            "N. gonorrhoeae", "N. meningitidis", "S. pneumoniae", "K. pneumoniae", "S. aureus", "P. aeruginosa"
        ]
    )

    _, xticks = plt.xticks()
    p.set_xticklabels(xticks, size = 14)
    p.set_yticklabels([round(f, 2) for f in p.get_yticks()], size = 14)
    plt.title("MLST prediction at 1000 reads [Methods]", size=16, fontweight="bold")
    plt.xlabel("")
    plt.ylabel("Mean proportion correctly predicted [n = 200, standard deviation]", size= 16)
    plt.legend(title='Species', fontsize='x-large', title_fontsize='16')
    plt.tight_layout()
    fig2.savefig(f"species_method_1000.png")

def _get_mlst_summary(eval_dir: Path, method: str = 'flye'):

    full_data = []
    for file in eval_dir.glob(f"*_{method}.txt"):
        species, _ = file.stem.split("_")
        ref = pandas.read_csv(eval_dir / f"{species}_truth.tsv", sep='\t', header=0)
        df = pandas.read_csv(file, sep='\t', header=None)
        names = [Path(name).stem.split("_")[0] for name in df.iloc[:, 0].tolist()]
        rep = [Path(name).stem.split("_")[1] for name in df.iloc[:, 0].tolist()]
        mlst = df.iloc[:, 2].tolist()
        mlst_sum = pandas.DataFrame({'name': names, 'replicate': rep, 'prediction': mlst})
        
        print(species, mlst_sum['prediction'].tolist())
        species_mlst = pandas.merge(mlst_sum, ref, left_on="name", right_on="accession")
        
        assert len(species_mlst) == len(mlst_sum)
        species_mlst['species'] = [SPECIES_DICT[species] for _ in species_mlst.iterrows()]
        full_data.append(species_mlst[['name', 'replicate', 'prediction', 'mlst', 'species']])

    df = pandas.concat(full_data)
    df['mlst'] = df['mlst'].astype(str)
    df['correct' ]= [
        1 if truth == prediction else 0 for truth, prediction in zip(df['mlst'], df['prediction'])
    ]

    print(f"Total correct {method}: {df['correct'].sum()}")

    summary_df = []
    for species, spec_data in df.groupby("species"):
        for rep, rep_data in spec_data.groupby("replicate"):
            if method == 'medaka':
                m = 'flye+medaka'
            else:
                m = 'flye'
            summary_df.append([species, rep, rep_data['correct'].sum() / len(rep_data), m])
    
    df = pandas.DataFrame(summary_df, columns=['species', 'rep', 'correct', 'method'])
    print(df)
    return df

def _get_krocus_summary(eval_dir: Path):

    krocus_data = []
    for f in eval_dir.glob("*_krocus.txt"):
        species, _ = f.stem.split("_")
        df = pandas.read_csv(f, sep=" ", header=None, names=['name', 'prediction'], na_values=["ND"])
        df['species'] = [SPECIES_DICT[species] for _ in df.iterrows()]
        df['method'] = ["krocus" for _ in df.iterrows()]
        krocus_data.append(df)

    df2 = pandas.concat(krocus_data)

    # Krocus did not return any predictions, set to 0
    df2 = df2.fillna(0)
    df2 = df2.rename(columns={'prediction': 'correct'})

    return df2


def _get_sketchy_summary(eval_dir: Path):

    full_data = []
    for file in list(eval_dir.glob("*_refdb.tsv")) + list(eval_dir.glob("*_simdb.tsv")):
        species, db = file.stem.split("_")
        df = pandas.read_csv(file, sep='\t')
        df['species'] = [species for _ in df.iterrows()]
        df['db'] = [db for _ in df.iterrows()]
        # offline for these plots
        df = df[df['mode'] == 'offline']
        full_data.append(df)

    df1 = pandas.concat(full_data)
    df1['correct'] = [
        1 if truth == prediction else 0 for truth, prediction in zip(df1['truth'], df1['prediction'])
    ]

    df_limit = df1[df1['reads'] == 1000]
    summary_1000 = []
    for species, spec_data in df_limit.groupby("species"):
        for db, db_data in spec_data.groupby("db"):
            for replicate, rep_data in db_data.groupby("replicate"):
                prop_correct = sum(rep_data['correct']) / len(rep_data)
                summary_1000.append([SPECIES_DICT[species], db, replicate, prop_correct])
    summary_1000 = pandas.DataFrame(summary_1000, columns=['species', 'db', 'replicate', 'correct'])

    return df1, summary_1000


@app.command()
def simsum(
    ref: Path,
    eval_dir: Path, 
    sim_sketch: Optional[bool] = typer.Option(None, "-s")
):

    ref_df = pandas.read_csv(ref, sep="\t", header=0)
    if sim_sketch:
        _name = 'simdb'
        data = []
        for file in eval_dir.glob("*.txt"):
            if "online" in file.name:
                name, reads, rep = file.name.strip(".online.txt").split("_")
                mode = 'online'
            else:
                name, reads, rep =  file.name.strip(".offline.txt").split("_")
                mode = 'offline'

            with file.open() as fin:
                _, mlst = fin.readline().strip().split(' ')
            
            lookup = ref_df.loc[ref_df['accession'] == name, 'mlst'].values[0]
            print(name, reads, mlst, lookup)
            data.append([name, int(reads), int(rep), int(mlst), int(lookup), mode])
        
        df = pandas.DataFrame(data, columns=['name', 'reads', 'replicate', 'prediction', 'truth', 'mode']).sort_values(['mode', 'replicate', 'reads'])
        df.to_csv(f"summary_{_name}.tsv", sep='\t', index=False)

    else:
        # ref sketch, offline only
        data = []
        _name = 'refdb'
        for file in eval_dir.glob("*.txt"):
            _, reads, rep = file.name.strip(".offline.txt").split("_")
            df = pandas.read_csv(file, sep=" ", header=None, names=['name', 'prediction']).sort_values('name')
            df['reads'] = [reads for _ in df.iterrows()]
            rdf = ref_df[ref_df['accession'].isin(df['name'])].sort_values('accession')
            assert df['name'].tolist() == rdf['accession'].tolist()
            df['truth'] = rdf['mlst'].tolist()
            df['mode'] = ['offline' for _ in df.iterrows()]
            df['replicate'] = [int(rep) for _ in df.iterrows()]
            data.append(df)
        df = pandas.concat(data).sort_values(['reads', 'name'])
        df.to_csv(f"summary_{_name}.tsv", sep='\t', index=False)

    results_data = []
    for mode, mdata in df.groupby('mode'):
        for read, read_data in mdata.groupby("reads", sort=True):
            read_truths = []
            for rep, rep_data in read_data.groupby("replicate"):
                truth_count = sum([1 if tup[0] == tup[1] else 0 for tup in zip(
                    rep_data['prediction'].tolist(), rep_data['truth'].tolist()
                )])
                read_truths.append(truth_count/len(rep_data))
            results_data.append([mode, int(read), mean(read_truths), stdev(read_truths)])
            print(mode, int(read), mean(read_truths), stdev(read_truths))

    rdf = pandas.DataFrame(results_data, columns=['mode', 'reads', 'mean', 'stdev']).sort_values(['mode', 'reads'])
    rdf['db'] = [_name for _ in rdf.iterrows()]
    rdf.to_csv(f"summary_results_{_name}.tsv", sep='\t', index=False)
    print(rdf)


@app.command()
def subset(
    ref: Path,
    dbs: Path,
    ext: Optional[str] = typer.Option(".txt", "-e")
):
    """
    Subset a genotype file based on reference sketch indices (sketchy info)
    """

    ref_df = pandas.read_csv(ref, sep="\t", header=0)

    dbidx = [
        f for f in dbs.glob(f"*{ext if ext else ''}") if f.is_file()
    ]

    for idx in dbidx:
        name = idx.stem
        print(name)
        df = pandas.read_csv(idx, sep=" ", header=None, names=['id', 'bp', 'unique'])
        out_df = ref_df[ref_df['accession'].isin(df['id'])]
        assert len(out_df) == len(df), 'length of data subset does not match genotype index length'
        out_df[['accession', 'mlst']].to_csv(f"{name}.tsv", sep='\t', index=False)


if __name__ == "__main__":
    app()