import click

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--metrics',
    '-m',
    type=Path,
    required=True,
    help='Path to directory with Fastq [required]',
)
@click.option(
    '--mpl_backend',
    type=str,
    default="",
    help='Matplotlib backend [default]'
)
@click.option(
    '--prefix',
    type=str,
    default="",
    help='Output plot prefix [none]'
)
def bootstrap_metrics(
    metrics, mpl_backend, prefix
):

    """ Barplot plot of Fastq read counts in a directory """

    sd = SketchyDiagnostics(outdir=None, mpl_backend=mpl_backend)
    sd.plot_bootstrap_metrics(metrics=metrics, prefix=prefix)

