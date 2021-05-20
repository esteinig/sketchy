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
    default="nfx-results",
    help='Nextflow summary output directory'
)
@click.option(
    '--exclude',
    '-e',
    required=False,
    type=str,
    default="",
    help='Run samples to exclude'
)
@click.option(
    '--force_db',
    '-f',
    required=False,
    type=str,
    default="",
    help='Force database [none]'
)
def match(
    directory, outdir, reference, exclude, force_db
):

    """ Match a reference genotype table to prediction outputs from Nextflow """

    to_exclude = [s.strip() for s in exclude.split(",")]

    sd = SketchyDiagnostics(outdir=outdir, verbose=True, mpl_backend=None)
    sd.match_reference(nextflow=directory, reference=reference, exclude_isolates=to_exclude, force_db=force_db)
