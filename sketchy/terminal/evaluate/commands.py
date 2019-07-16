import click
import pandas
import shutil

from sketchy.evaluation import SampleEvaluator
from pathlib import Path


@click.command()
@click.option(
    '--indir', '-d', default=None, required=True,  type=Path,
    help='Input directory from sketchy predict --keep.'
)
@click.option(
    '--truth', '-t', default=None, required=False,  type=Path,
    help='Path to JSON file with truth keys: lineage, susceptibility, genotype.'
)
@click.option(
    '--outdir', '-o', default='sample_evaluation', type=Path,
    help='Output directory for evaluation data and plots.'
)
@click.option(
    '--limit', '-l', default=None,  type=int,
    help='Evaluate up to and including this number of reads.'
)
def evaluate(indir, truth, outdir, limit):

    """ Evaluate a sample for detection boundaries """

    SampleEvaluator(indir, outdir, limit=limit)
