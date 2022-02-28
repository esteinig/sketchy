import pandas
from pathlib import Path


def compute_metrics(
    reference: Path,
    predictions: Path = None,
    sketch_genotypes: Path = None,
    rase_prediction: Path = None,
    output: Path = Path("metrics.tsv")
):
    """
    Computes summary metrics for predictions against a structured reference

    :params reference: tab-delineated file with headers containing reference classifications
    :params predictions: tab-delineated file with headers containing prediction classfications
        must be the same order as `reference` table
    :params sketch_genotypes"
    """

    if predictions is not None:
        ref_df, prediction_df = process_data(reference, predictions, sketch_genotypes)


def process_data(reference_genotypes: Path, predictions: Path, sketch_genotypes: Path):
    """
    Process reference and prediction tables, ensure that data is clean
    """

    ref_df = read_reference_genotypes(file=reference_genotypes)
    prediction_df = read_predictions(file=predictions, genotype_file=sketch_genotypes)

    # Retain only predictions that are also in the reference genotype data
    prediction_df = prediction_df[prediction_df['name'].isin(ref_df['name'])]

    # Retain only references that are also in the predictions
    ref_df = ref_df[ref_df['name'].isin(prediction_df['name'])]

    # Sort the dataframes in the same order by name
    ref_df = ref_df.sort_values('name')
    prediction_df = prediction_df.sort_values('name')

    assert ref_df['name'].tolist() == prediction_df['name'].tolist()

    # For the S. aureus data predictions can be PVL* (which is PVL negative, missing one gene)
    # and resistances can include 'r' instead of 'R' - harmonize these by replacing with
    # appropriate value (r -> R, PVL* -> PVL-)

    ref_df = ref_df.replace('r', 'R')
    ref_df = ref_df.replace('PVL*', 'PVL-')

    prediction_df = prediction_df.replace('r', 'R')
    prediction_df = prediction_df.replace('PVL*', 'PVL-')

    # Also remove leading whitespace from SCCmec types
    ref_df['scc'] = ref_df.scc.str.strip()
    prediction_df['scc'] = prediction_df.scc.str.strip()

    # Drop mecA gene (assembly based) in favour of Mykrobe
    # methicillin typing from reads (same genotype)
    ref_df = ref_df.drop(columns='meca')
    prediction_df = prediction_df.drop(columns='meca')
    ref_df = ref_df.reset_index(drop=True)
    prediction_df = prediction_df.reset_index(drop=True)

    return ref_df, prediction_df


def read_reference_genotypes(file: Path):
    """
    Read the reference genotype file, with headers
    """
    return pandas.read_csv(file, sep="\t", header=0)


def read_predictions(file: Path, genotype_file: Path):
    sketch_genotypes = pandas.read_csv(genotype_file, sep='\t')
    genotype_columns = sketch_genotypes.columns.tolist()[1:]  # removes id column

    cols = ["name", "reads", "sketch_id", "shared_hashes"] + genotype_columns

    return pandas.read_csv(
        file, sep='\s+', names=cols, converters={"name": lambda x: Path(x).stem}
    )  # whitespace or tab delimiters