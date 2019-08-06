import click

from pathlib import Path
from sketchy.minhash import MashSketch


@click.command()
@click.option(
    '--data', '-d', help='Input data file to assign lineages and genotypes.'
)
@click.option(
    '--outdir', '-o', default=Path().cwd(), type=Path,
    help='Output directory for sketch files.'
)
@click.option(
    '--copy', '-c', is_flag=True,
    help='Copy assembly files.'
)
@click.option(
    '--delimiter', default='\t',
    help='Delimiter for data file.'
)
@click.option(
    '--uuid', is_flag=True,
    help='Use a UUID instead of the data index to anonymise data collection.'
)
def sk_link(data, outdir, copy, delimiter, uuid):
    """ Symlink a set of .fasta files from file """

    sketch = MashSketch()

    sketch.read_data(
        fpath=data, sep=delimiter, index=0
    )

    _ = sketch.link(
        fdir=outdir,
        rename=None,
        symlink=not copy,
        progbar=True,
        uuid=uuid
    )



