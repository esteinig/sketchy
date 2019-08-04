import click
import pandas

from sketchy.evaluation import SampleEvaluator
from pathlib import Path
from matplotlib import pyplot as plt

@click.command()
@click.option(
    '--indir', '-i', default=None, required=True,  type=Path,
    help='Input directory from sketchy predict --keep.'
)
@click.option(
    '--data', '-d',  type=str, default=None, required=True,
    help='MASH sketch data; or a template, one of: kleb, mrsa, tb'
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
    '--ranks', default=100,  type=int,
    help='Secondary color for hitmap (lineage correct only).'
)
@click.option(
    '--top', default=50,  type=int,
    help='Collect the top ranked genome hits by sum of shared hashes to plot.'
)
@click.option(
    '--top_lineages', default=5,  type=int,
    help='Collect the top ranked lineages aggregated by sum of sums of '
         'shared hashes for plotting.'
)
@click.option(
    '--multi', '-m', is_flag=True,
    help='Output parsed is raw MASH output from prediction with --ncpu > 1'
)
@click.option(
    '--sketchy',  default=Path.home() / '.sketchy', type=Path,
    help='Path to Sketchy home directory [ ~/.sketchy/ ]'
)
def evaluate(
    indir, lineage, resistance, genotype, outdir, limit, multi, data,
    color, primary, secondary, ranks, top, top_lineages, sketchy
):

    """ Evaluate a sample for detection boundaries """
    if data in ('kleb', 'mrsa', 'tb'):
        sketch_data = pandas.read_csv(
            sketchy / 'data' / f'{data}.data.tsv', sep='\t', index_col=0
        )
    else:
        sketch_data = pandas.read_csv(
            Path(data), sep='\t', index_col=0
        )

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
        sequential=not multi,
        sketch_data=sketch_data
    )

    validate = [lineage, resistance, genotype]
    if validate.count(None) == len(validate):
        fig, (ax1, ax2) = plt.subplots(
            nrows=1, ncols=2, figsize=(21.0, 7.0)
        )
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle(f'{indir.name} - {limit}')

        se.create_lineage_hitmap(top=top_lineages, ax=ax1)
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

        se.create_validation_hitmap(ranks=ranks, ax=ax1)
        se.create_race_plot(ax=ax2)
        se.create_concordance_plot(ax=ax3)

        fig.savefig(
            outdir / 'validation_plots.pdf',
        )