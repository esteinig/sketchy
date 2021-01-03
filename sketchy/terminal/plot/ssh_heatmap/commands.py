import click

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--db',
    '-d',
    type=Path,
    required=True,
    help='Path to database directory used in the run; full path required [required]',
)
@click.option(
    '--ssh',
    '-s',
    type=Path,
    required=True,
    help='Path to raw sums of shared hashes data file from stream [required]',
)
@click.option(
    '--outdir',
    '-o',
    type=Path,
    default=Path("diagnostics"),
    help='Output directory for diagnostic files and plots [diagnostics]'
)
@click.option(
    '--plot',
    '-p',
    type=Path,
    default=Path("ssh.png"),
    help='Plot file, extension specifies format [ssh.png]'
)
@click.option(
    '--color',
    '-c',
    type=str,
    default='YlGnBu',
    help='Color palette for output plots [YlGnBu]'
)
@click.option(
    '--max_ranks',
    '-m',
    type=int,
    default=5,
    help='Number shared hashes ranks to color; others are grayed out [5]'
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
def ssh_heatmap(db, ssh, plot, max_ranks, color, outdir, mpl_backend, verbose):

    """ Diagnostic heatmap for ranked sum of shared hashes """

    sd = SketchyDiagnostics(outdir=outdir, verbose=verbose, mpl_backend=mpl_backend)
    sd.plot_ssh_diagnostics(db=db, ssh_file=ssh, plot_file=plot, max_ranks=max_ranks, color=color)



