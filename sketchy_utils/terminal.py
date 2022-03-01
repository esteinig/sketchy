import typer
from pathlib import Path
from .metrics import compute_metrics

app = typer.Typer(add_completion=False)


@app.command()
def metrics(
    reference: Path = typer.Argument(
        None, help="Reference genotype file with `name` column and genotypes in the same format as predictions"
    ),
    prediction: Path = typer.Argument(
        None, help="Prediction genotype file with `name` column and outputs from `sketchy predict`"
    ),
    label_config: Path = typer.Option(
        None, help="Label configuration file (JSON) to specify whether a prediction feature (column) is binary"
                   "or multi-label, and to specify whether non standard binary labels should be used. "
                   "Prediction column names as keys, and dictionaries as values: `binary`: bool, "
                   "`labels`: list (empty or binary labels). Labels only need to be specified for binary columns"
                   "and an empty list will use the default binary option as replacement."
    ),
    default_binary: str = typer.Option(
        "R,S", help="Default binary labels, in a comma-separated string [e.g. 'R,S'] used for all binary features"
                    "where labels are not specified in the configuration, or for all features not present in "
                    "the configuration (default feature is binary)"
    ),
    output: Path = typer.Option(
        Path("metrics.tsv"), help="Output summary file containing computed feature metrics"
    ),
    verbose: bool = typer.Option(
        False, help="Verbose outputs to make sure features are correctly assigned binary status and labels"
    )
):

    compute_metrics(
        reference=reference,
        prediction=prediction,
        label_config=label_config,
        default_binary=default_binary.split(","),
        output=output,
        verbose=verbose
    )
