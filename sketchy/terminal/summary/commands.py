import click

from pathlib import Path
from sketchy.minhash import MashScore


@click.command()
@click.option(
    '--indir', '-i', type=Path, help='Temporary directory with read hashes (--ncpu > 0).'
)
@click.option(
    '--outdir', '-o', type=Path, help='Output directory for read-wise sum of shared hashes.'
)
def summary(indir, outdir):
    """ Sum a temporary read hashes into the sum of shared hashes """

    mash_score = MashScore()
    mash_score.sum_read_hashes(
        indir=indir, outdir=outdir
    )
