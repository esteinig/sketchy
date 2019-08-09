import click
import pandas

from sketchy.evaluation import SampleEvaluator
from matplotlib import pyplot as plt
from pathlib import Path


@click.command()
@click.option(
    '--data', '-d', help='Sketchy data file from predict output', type=Path, required=True,
)
@click.option(
    '--top', '-t', default=5,  type=int,
    help='Select most frequent lineages to color and plot'
)
@click.option(
    '--genotype', '-g', is_flag=True,
    help='Enable genotype sub-plots'
)
@click.option(
    '--resistance', '-r', is_flag=True,
    help='Enable resistance profile sub-plots'
)
@click.option(
    '--breakpoints', '-b', is_flag=True,
    help='Enable breakpoint computation for first and '
         'stable detection of most common trait.'
)
@click.option(
    '--stable', '-s', default=300, type=int,
    help='Define block size of consecutive prediction'
         'across reads to define a stable breakpoint.'
)
@click.option(
    '--color', '-c', default='BuGn',
    help='Color palette, Brewer colors for Seaborn'
)
@click.option(
    '--format', '-f', default='png',  type=str,
    help='Output format for plot file'
)
@click.option(
    '--prefix', '-p', default='sketchy',  type=str,
    help='Prefix for plot file'
)
@click.option(
    '--limit', '-l', default=None,  type=int,
    help='Limit the number of reads along x-axis in plots.'
)
def plot(data, top, genotype, resistance, limit, format, prefix, color, breakpoints, stable):
    """ Generate lineage hitmap and sum plots from predictions. """

    df = pandas.read_csv(data, sep='\t', index_col=0)

    split_colors = color.split(',')

    if len(split_colors) == 3:
        color, gcolor, rcolor = split_colors
    elif len(split_colors) == 2:
        color, gcolor = split_colors
        rcolor = None
    else:
        color, gcolor, rcolor = color, color, color

    nreads = len(
        df['read'].unique()
    )

    se = SampleEvaluator(
        limit=limit,
        plot_data=df,
        top=len(
            df['rank'].unique()
        ),
    )

    if limit and limit > nreads:
        se.logger.info(
            f'Selected limit {limit} is larger than number of reads {nreads}'
        )
        exit(1)

    if genotype and not resistance:
        nrows, ncols = 2, 2
    elif resistance and not genotype:
        nrows, ncols = 2, 2
    elif resistance and genotype:
        nrows, ncols = 3, 2
    else:
        nrows, ncols = 1, 2

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(14.0, nrows*4.5)
    )

    fig.subplots_adjust(hspace=0.8)
    fig.suptitle('')

    se.logger.info(f'Generate hitmap [ {top}, {se.limit} ] for lineage ...')
    se.logger.info(
        f'Generate lineplot for sum of sums of shared hashes by lineage ...'
    )

    se.create_hitmap(top=top, ax=axes[0, 0], color=color)
    se.create_lineplot(top=top, ax=axes[0, 1], color=color)

    if genotype and resistance:
        by = 'genotype and susceptibility'
    else:
        by = 'genotype' if genotype and not resistance else 'susceptibility'

    se.logger.info(f'Generate hitmap [ {top}, {se.limit} ] for {by} ...')
    se.logger.info(
        f'Generate lineplot for sum of sums of shared hashes by {by} ...'
    )

    if genotype and not resistance:
        se.create_hitmap(
            top=top, data='genotype', ax=axes[1, 0], color=gcolor
        )
        se.create_lineplot(
            top=top, data='genotype', ax=axes[1, 1], color=gcolor
        )
    elif resistance and not genotype:

        se.create_hitmap(
            top=top, data='susceptibility', ax=axes[1, 0], color=rcolor
        )
        se.create_lineplot(
            top=top, data='susceptibility', ax=axes[1, 1], color=rcolor
        )
    elif resistance and genotype:
        se.create_hitmap(
            top=top, data='genotype', ax=axes[1, 0], color=gcolor
        )
        se.create_lineplot(
            top=top, data='genotype', ax=axes[1, 1], color=gcolor
        )
        se.create_hitmap(
            top=top, data='susceptibility', ax=axes[2, 0], color=rcolor
        )
        se.create_lineplot(
            top=top, data='susceptibility', ax=axes[2, 1], color=rcolor
        )

    if breakpoints:

        if stable > se.limit:
            se.logger.info(
                f'Breakpoints could not be computed: block detection {stable} '
                f'(--stable) must be smaller than last evaluated read @ {se.limit}'
            )
            se.logger.info(
                f'Exiting ...'
            )  # Nextflow break
            exit(1)

        b1 = se.find_breakpoints(
            top=top, data='lineage', block_size=stable
        )

        b2, b3 = None, None
        if genotype:
            b2 = se.find_breakpoints(
                top=top, data='genotype', block_size=stable
            )
        if resistance:
            b3 = se.find_breakpoints(
                top=top, data='susceptibility', block_size=stable
            )

        bp = pandas.DataFrame(
            data={
                'lineage': b1,
                'genotype': b2,
                'susceptibility': b3
            },
            index=['first', 'stable', 'prediction']
        )

        bp.to_csv(f'{prefix}.bp.tsv', sep='\t')

    plt.tight_layout()

    fig.savefig(
        f'{prefix}.{format}'
    )
