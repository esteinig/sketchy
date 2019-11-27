import click

from pathlib import Path
from sketchy.tests.minhash import MashSketch


@click.command()
@click.option(
    '--fasta', '-f', type=Path,
    help='Fasta (.fasta) input directory to sketch genomes.'
)
@click.option(
    '--prefix', '-', default='sketch',
    help='Prefix for sketch output.'
)
@click.option(
    '--kmer', '-k', default=15,
    help='K-mer length in MASH.'
)
@click.option(
    '--size', '-s', default=1000,
    help='Sketch size in MASH.'
)
@click.option(
    '--glob', '-g', default="*",
    help='Glob for Fasta files in directory --fasta.'
)
def sk_sketch(fasta, kmer, prefix, glob, size):
    """Create a MinHash sketch, wraps MASH"""

    sketch = MashSketch()

    sketch.sketch(
        name=Path(f'{prefix}_{kmer}_{size}'),
        fdir=fasta,
        k=kmer,
        size=size,
        glob=glob
    )



