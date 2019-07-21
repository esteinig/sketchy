import click

from sketchy.evaluation import SampleEvaluator
from pathlib import Path


@click.command()
@click.option(
    '--indir', '-d', default=None, required=True,  type=Path,
    help='Input directory from sketchy predict --keep.'
)
@click.option(
    '--outdir', '-o', default='sample_evaluation', type=Path,
    help='Output directory for evaluation data and plots.'
)
@click.option(
    '--limit', '-l', default=1000,  type=int,
    help='Evaluate up to and including this number of reads.'
)
@click.option(
    '--color', '-c', default=None,  type=str,
    help='Color of heatmap output: red,orange,green,blue'
)
@click.option(
    '--lineage', default='9',  type=str,
    help='True lineage to evaluate on.'
)
@click.option(
    '--resistance', default='SRSSSSSSRSSS',  type=str,
    help='True resistance profile to evaluate on.'
)
@click.option(
    '--genotype', default='',  type=str,
    help='True genotype to evaluate on.'
)
@click.option(
    '--primary', default="#88419d",  type=str,
    help='Primary color for hitmap (joint lineage, genotype, susceptibility).'
)
@click.option(
    '--secondary', default="#8c96c6",  type=str,
    help='Secondary color for hitmap (lineage correct only).'
)
def evaluate(indir, lineage, resistance, genotype, outdir, limit, color, primary, secondary):

    """ Evaluate a sample for detection boundaries """

    se = SampleEvaluator(
        indir, outdir,
        limit=limit,
        palette=color,
        true_lineage=lineage,
        true_resistance=resistance,
        true_genotype=genotype,
        primary_color=primary,
        secondary_color=secondary,
    )

    se.create_timeline_hitmap()
    se.create_race_plot()
    se.create_concordance_plot()
