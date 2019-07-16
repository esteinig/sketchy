import click
import pandas
import shutil

from sketchy.evaluation import BootstrapEvaluator
from pathlib import Path


@click.command()
@click.option(
    '--file', '-f', default=Path().cwd(), help='Bootstrap prediction files of all bootstraps.'
)
@click.option(
    '--confidence', '-c', default=0.95, help='Confidence interval for bootstraps.'
)
def pplot(file, confidence):
    """ Plot a collection of bootstrap predictions. """

    ev = Evaluator(
        outdir=None
    )

    ev.plot_bootstraps(
        bootstrap_data=Path(file), confidence=confidence, display=False
    )


