import click

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--directory',
    '-d',
    type=Path,
    required=True,
    help='Path to directory with Fastq [required]',
)
@click.option(
    '--outdir',
    '-o',
    type=Path,
    required=False,
    default=Path.cwd(),
    help='Path to output directory for plot [cwd]',
)
@click.option(
    '--mpl_backend',
    type=str,
    default="",
    help='Matplotlib backend [default]'
)
@click.option(
    '--extension',
    type=str,
    default=".fastq",
    help='Fastq extension [.fastq]'
)
@click.option(
    '--prefix',
    type=str,
    default="",
    help='Output plot prefix [none]'
)
def barcode_counts(
    directory, mpl_backend, outdir, extension, prefix
):

    """ Barplot plot of Fastq read counts in a directory """

    sd = SketchyDiagnostics(outdir=outdir, mpl_backend=mpl_backend)
    sd.plot_barcode_barplot(directory=directory, ext=extension, prefix=prefix)


