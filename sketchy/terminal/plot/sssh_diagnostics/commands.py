import click

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--sssh',
    '-s',
    type=Path,
    required=True,
    help='Path to sum of ranked sums shared hashes data file from evaluation',
)
@click.option(
    '--outdir',
    '-o',
    type=Path,
    default=Path("diagnostics"),
    help='Output directory for diagnostic files and plots [diagnostics]'
)
@click.option(
    '--plot_file',
    '-p',
    type=Path,
    default=Path("ssh.png"),
    help='Plot file, extension specifies format [diagnostics.png]'
)
@click.option(
    '--stable',
    '-b',
    type=int,
    default=100,
    help='Stability parameter passed to: sketchy stream to compute stable breakpoint for each feature [none]'
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
    help='In the plots, number of feature values / prediction ranks; reduces cluttering '
         'if many alternatives genotypes called [5]'
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
def sssh_diagnostics(sssh, plot_file, stable, max_ranks, color, outdir, mpl_backend, verbose):

    """ Sum of ranked sums of shared hashes diagnostics (main) """

    sd = SketchyDiagnostics(outdir=outdir, verbose=verbose, mpl_backend=mpl_backend)

    sssh_data = sd.process_sssh(sssh_file=sssh, stable=stable, max_ranks=max_ranks, mode="last")
    sd.plot_sssh_diagostics(sssh_data=sssh_data, plot_file=plot_file, color=color)

