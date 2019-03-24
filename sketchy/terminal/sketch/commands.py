import click

from pathlib import Path
from sketchy.sketchy import Sketchy


@click.command()
@click.option(
    '--data', '-d', help='Input data file to assign lineages and genotypes.'
)
@click.option(
    '--kmer', '-k', default=15,
    help='K-mer length to create sketch with in MASH.'
)
@click.option(
    '--size', '-s', default=1000,
    help='Sketch size in MASH.'
)
@click.option(
    '--outdir', '-o', default=Path().cwd() / 'sketchy_sketch',
    help='output directory for sketch files.'
)
@click.option(
    '--prefix', '-p', default='sketchy',
    help='Prefix for data sketch with Sketchy.'
)
@click.option(
    '--copy', '-c', is_flag=True,
    help='Copy assembly files.'
)
@click.option(
    '--delimiter', default='\t',
    help='Delimiter for data file.'
)
def sketch(data, kmer, outdir, copy, prefix, delimiter, size):
    """Create a MinHash sketch for matching with MASH"""

    sketchy = Sketchy()
    sketchy.create_mash_sketch(
        data=Path(data).resolve(),
        outdir=Path(outdir).resolve(),
        kmer_length=kmer,
        sketch_size=size,
        file_glob='*',
        file_copy=copy,
        prefix=prefix,
        sep=delimiter
    )


