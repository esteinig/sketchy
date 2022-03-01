import pandas
import json
from numpy import nan
from typing import List
from pathlib import Path
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import confusion_matrix

SKIP_COLUMNS = ["name", "reads", "sketch_id", "shared_hashes"]


def compute_metrics(
    reference: Path,
    prediction: Path = None,
    label_config: Path = None,
    default_binary: List = ("R", "S"),
    output: Path = Path("metrics.tsv"),
    verbose: bool = True
):
    """
    Computes summary metrics for predictions against a structured reference

    :params reference: a tab-delineated file with headers containing reference classifications
        * must contain a column `name` for file name matching with prediction file

    :params predictions: a tab-delineated file with headers containing prediction classfications
        * must contain a column 'name' for file name same as in reference file
        * this can be the output from the batch prediction pipeline

    :params label_config: a JSON file in the following structure:

        {
            "<column_name>": {
                "binary": true/false,
                "labels": [<label1>, <label2>] or null
        }

        If <binary> is true, then a list of two labels can be provided, which correspond to the
        postive and negative values for the confusion matrix.

        {
            "mlst": {
                "binary": false,
                "labels": null
            },
            "pvl": {
                "binary": true,
                "labels": ["PVL+", "PVL-"]
            },
            "mrsa": {
                "binary": true,
                "labels": ["MRSA", "MSSA"]
            }
        }

    """

    reference_df, prediction_df = process_data(reference=reference, prediction=prediction, verbose=verbose)

    if label_config is not None:
        config = read_label_config(file=label_config, default_binary=default_binary)
    else:
        config = {}

    metrics = get_metrics(
        reference=reference_df, prediction=prediction_df, config=config, default_binary=default_binary, verbose=verbose
    )

    metrics.to_csv(output, sep='\t', index=False)
    print(metrics)


def process_data(reference: Path, prediction: Path, verbose: bool):
    """
    Process reference and prediction tables, ensure that data is clean
    """

    ref = pandas.read_csv(reference, sep="\t", header=0)
    pred = pandas.read_csv(prediction, sep="\t", header=0)

    # Retain only predictions that are also in the reference genotype data
    prediction_df_clean = pred[pred['name'].isin(ref['name'])]

    # Retain only references that are also in the predictions
    ref_df_clean = ref[ref['name'].isin(pred['name'])]

    removed_prediction = len(pred) - len(prediction_df_clean)
    removed_ref = len(ref) - len(ref_df_clean)

    if verbose:
        print(f"Removed {removed_prediction} entries from prediction data not present in reference data")
        print(f"Removed {removed_ref} entries from reference data not present in prediction data")

    # Sort the dataframes in the same order by name
    ref_df = ref_df_clean.sort_values('name')
    prediction_df = prediction_df_clean.sort_values('name')

    assert ref_df['name'].tolist() == prediction_df['name'].tolist()  # names should be unique, therefore sortable

    # Make sure all columns in prediction are present in reference:

    ref_columns = ref_df.columns.tolist()
    for column in prediction_df:
        if column in SKIP_COLUMNS:
            continue
        else:
            if column not in ref_columns:
                raise ValueError(f"Column `{column}` not in reference data")

    # Following should be done in the reference genotype files later on!

    # For the S. aureus data predictions can be PVL* (which is PVL negative, missing one gene)
    # and resistances can include 'r' instead of 'R' - harmonize these by replacing with
    # appropriate value (r -> R, PVL* -> PVL-)

    ref_df = ref_df.replace('r', 'R')
    prediction_df = prediction_df.replace('r', 'R')

    if "pvl" in ref_df.columns and "pvl" in prediction_df.columns:
        ref_df = ref_df.replace('PVL*', 'PVL-')
        prediction_df = prediction_df.replace('PVL*', 'PVL-')

    if "scc" in ref_df.columns and "scc" in prediction_df.columns:
        # Also remove leading whitespace from SCCmec types
        ref_df['scc'] = ref_df.scc.str.strip()
        prediction_df['scc'] = prediction_df.scc.str.strip()

    if "meca" in ref_df.columns and "meca" in prediction_df.columns:
        # Drop mecA gene (assembly based) in favour of Mykrobe
        # methicillin typing from reads (same genotype)
        ref_df = ref_df.drop(columns='meca')
        prediction_df = prediction_df.drop(columns='meca')
        ref_df = ref_df.reset_index(drop=True)
        prediction_df = prediction_df.reset_index(drop=True)

    return ref_df, prediction_df


def read_label_config(file: Path, default_binary: list) -> dict:
    """
    Reads and validates the label configuration JSON
    """

    with file.open() as json_file:
        label_config = json.load(json_file)

    config = {}
    for column, config in label_config.items():
        if "binary" not in config.keys():
            raise ValueError(f"Could not find `binary` key in config for column: {column}")
        if "labels" not in config.keys():
            raise ValueError(f"Could not find `labels` key in config for column: {column}")

        binary: bool = config["binary"]

        if not isinstance(binary, bool):
            raise ValueError(f"Binary value ({column}) is not a boolean")

        config[column] = {"binary": binary}

        if binary:
            labels: list = config["labels"]
            if labels is not None or not isinstance(labels, list):
                raise ValueError(f"Label value ({column}) is not a list or None")

            if labels is None:
                print(f"Setting default binary labels ({column}): {default_binary}")
                labels = default_binary.copy()
            else:
                if len(labels) != 2:
                    raise ValueError(f"Binary labels requires precisely two values ({column})")

            config[column]["labels"] = labels

    return config


def get_metrics(
    reference: pandas.DataFrame, prediction: pandas.DataFrame, config: dict, default_binary: List, verbose: bool
) -> pandas.DataFrame:

    metrics = []
    for column in prediction.columns:
        # Skip irrelevant columns
        if column in SKIP_COLUMNS:
            continue
        else:

            try:
                binary = config[column]["binary"]
                if verbose:
                    print(f"Column ({column}) is binary: {binary}")
            except KeyError:
                if verbose:
                    print(f"Column ({column}) is binary (default)")
                binary = True

            ref_vec = reference[column].tolist()  # sorted in process_data
            pred_vec = prediction[column].tolist()

            if binary:
                try:
                    labels = config[column]["labels"]
                    if verbose:
                        print(f"Column ({column}) labels: {', '.join(default_binary)}")
                except KeyError:
                    if verbose:
                        print(f"Column ({column}) labels: {', '.join(default_binary)} (default)")
                    labels = default_binary.copy()

                tp, fp, tn, fn, acc, tpr, tnr, ppv, npv = compute_binary_metrics(ref_vec, pred_vec, labels)
            else:
                tp, fp, tn, fn, tnr = nan, nan, nan, nan, nan
                acc = accuracy_score(ref_vec, pred_vec)
                tpr = recall_score(ref_vec, pred_vec, average='weighted')
                ppv = precision_score(ref_vec, pred_vec, average='weighted')

            metrics.append([column, binary, tp, tn, fp, fn, acc, ppv, tpr, tnr])

    return pandas.DataFrame(
        metrics,
        columns=[
            "feature", "binary", "true_positives", "true_negatives",
            "false_positives", "false_negatives", "accuracy", "precision",
            "recall", "specificity"
        ]
    )


def compute_binary_metrics(ref: list, pred: list, labels: list):
    try:
        cm = confusion_matrix(ref, pred, labels=labels)
    except ValueError:
        print(ref, pred, labels)
        raise

    tp = int(cm[0][0])
    fn = int(cm[0][1])
    fp = int(cm[1][0])
    tn = int(cm[1][1])

    # In all cases, if either the numerator or
    # denominator is zero, the metric is undefined
    tpr = _metric(tp, (tp + fn))
    tnr = _metric(tn, (tn + fp))
    ppv = _metric(tp, (tp + fp))
    npv = _metric(tn, (tn + fn))
    acc = _metric((tp + tn), (tp + fp + fn + tn))

    return tp, fp, tn, fn, acc, tpr, tnr, ppv, npv


def _metric(numerator: int, denominator: int):

    if numerator == 0:
        return nan
    if denominator == 0:
        return nan

    return numerator/denominator
