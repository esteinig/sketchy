import click

from sketchy.evaluation import SampleEvaluator
from pathlib import Path
from matplotlib import pyplot as plt

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
    '--lineage', default=None,  type=str or None,
    help='True lineage to evaluate on.'
)
@click.option(
    '--resistance', default=None,  type=str or None,
    help='True resistance profile to evaluate on.'
)
@click.option(
    '--genotype', default=None,  type=str or None,
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
@click.option(
    '--show_ranks', '--ranks', default=100,  type=int,
    help='Secondary color for hitmap (lineage correct only).'
)
@click.option(
    '--top', default=50,  type=int,
    help='Collect the top ranked genome hits by sum of shared hashes to plot.'
)
@click.option(
    '--top_lineages', default=10,  type=int,
    help='Collect the top ranked lineages aggregated by mean sum of '
         'shared hashes for plotting.'
)
def evaluate(indir, lineage, resistance, genotype, outdir, limit, color, primary, secondary, show_ranks, top, top_lineages):

    """ Evaluate a sample for detection boundaries """

    se = SampleEvaluator(
        indir, outdir,
        limit=limit,
        palette=color,
        top=top,
        true_lineage=lineage,
        true_resistance=resistance,
        true_genotype=genotype,
        primary_color=primary,
        secondary_color=secondary,
    )

    validate = [lineage, resistance, genotype]
    if validate.count(None) == len(validate):
        print('Computing plots...')
        fig, (ax1, ax2) = plt.subplots(
            nrows=1, ncols=2, figsize=(21.0, 7.0)
        )
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle(f'{indir.name} - {limit}')

        se.create_lineage_hitmap(top=top_lineages, ranks=show_ranks, ax=ax1)
        se.create_lineage_plot(top=top_lineages, ax=ax2)

        plt.tight_layout()

        fig.savefig(
            outdir / 'lineage_plots.pdf',
        )

    else:
        fig, (ax1, ax2, ax3) = plt.subplots(
            nrows=1, ncols=3, figsize=(21.0, 7.0)
        )
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle(f'{indir.name} - {limit}')

        se.create_timeline_hitmap(ranks=show_ranks, ax=ax1)
        se.create_race_plot(ax=ax2)
        se.create_concordance_plot(ax=ax3)

        fig.savefig(
            outdir / 'validation_plots.pdf',
        )