import click
from numpy import where

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--results',
    '-n',
    type=Path,
    required=True,
    help='Path to collected results directory of Nextflow run [required]',
)
@click.option(
    '--reference',
    '-r',
    type=Path,
    required=False,
    help='Path to reference genotype file for heatmap colors [optional]',
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
def genotype_heatmap(results, outdir, reference, plot, color, subset_column, subset_values, mpl_backend, verbose):

    """ Comparison of genotype predictions from Nextflow """

    sd = SketchyDiagnostics(outdir=outdir, verbose=verbose, mpl_backend=mpl_backend)

    if reference:
        sd.match_reference(nextflow=results, reference=reference)

    sd.plot_genotype_heatmap(
        nextflow=results, subset_column=subset_column, subset_values=subset_values
    )

