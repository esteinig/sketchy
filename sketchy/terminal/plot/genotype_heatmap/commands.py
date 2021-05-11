import click
from numpy import where

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--directory',
    '-d',
    type=Path,
    required=True,
    help='Path to collected results directory of Nextflow run [required]',
)
@click.option(
    '--outdir',
    '-o',
    type=Path,
    default=Path("nfx-heatmaps"),
    help='Output directory for plots [nfx-heatmaps]'
)
@click.option(
    '--plot',
    '-p',
    type=Path,
    default=Path("diagnostics.png"),
    help='Plot file, extension specifies format [diagnostics.png]'
)
@click.option(
    '--color',
    '-c',
    type=str,
    default='YlGnBu',
    help='Color palette for output plots [YlGnBu]'
)
@click.option(
    '--scale',
    '-s',
    type=float,
    default=1.0,
    help='Scale factor of plot size [1.0]'
)
@click.option(
    '--subset_column',
    '-sc',
    type=str,
    default=None,
    help='Column name to subset data [none]'
)
@click.option(
    '--subset_values',
    '-sv',
    type=str,
    default='None',
    help='Comma delimited string of values in subset column to select [none]'
)
@click.option(
    '--reverse_subset',
    '-rs',
    is_flag=True,
    help='Reverse the subset (selection NOT in column)'
)
@click.option(
    '--exclude_isolates',
    '-ei',
    type=str,
    default=None,
    help='Exclude isolates from heatmap output [none]'
)
@click.option(
    '--exclude_genotypes',
    '-eg',
    type=str,
    default=None,
    help='Exclude genotypes from heatmap output [none]'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose logging output [false]'
)
@click.option(
    '--mpl_backend',
    type=str,
    default="",
    help='Matplotlib backend [default]'
)
def genotype_heatmap(
    directory, outdir, scale, reverse_subset,
    subset_column, subset_values, mpl_backend, verbose, exclude_isolates, exclude_genotypes
):

    """ Comparison of genotype predictions from Nextflow """

    sd = SketchyDiagnostics(outdir=outdir, verbose=verbose, mpl_backend=mpl_backend)

    if exclude_genotypes:
        exclude_genotypes = [i.strip() for i in exclude_genotypes.split(',')]
    if exclude_isolates:
        exclude_isolates = [i.strip() for i in exclude_isolates.split(',')]

    sd.plot_genotype_heatmap(
        nextflow=directory, subset_column=subset_column, subset_values=subset_values, reverse_subset=reverse_subset,
        exclude_genotypes=exclude_genotypes, exclude_isolates=exclude_isolates, scale=scale
    )

