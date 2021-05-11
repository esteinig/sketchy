import click

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--directory',
    '-d',
    required=True,
    type=Path,
    help='Path to output directory from collected Nextflow'
)
@click.option(
    '--reference',
    '-r',
    type=Path,
    required=True,
    help='Path to reference genotype file for heatmap colors [optional]',
)
@click.option(
    '--outdir',
    '-o',
    required=False,
    type=Path,
    default="nxf-results",
    help='Nextflow summary output directory'
)
def match(
    directory, outdir, reference
):

    """ Match a reference genotype table to prediction outputs from Nextflow """

    sd = SketchyDiagnostics(outdir=outdir, verbose=True, mpl_backend=None)
    sd.match_reference(nextflow=directory, reference=reference)
