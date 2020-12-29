import click
import pandas
import tempfile
import numpy as np

from pathlib import Path
from sketchy.sketchy import LineageIndex
from sketchy.utils import run_cmd


@click.command()
@click.option(
    '--fasta', '-f', default=None, type=Path, required=True,
    help='Path to list containing path to genomes per line [required]'
)
@click.option(
    '--index', '-i', default=None, type=Path, required=True,
    help='Path to db_lineage index file [required]'
)
@click.option(
    '--output', '-o', default="db_lineage.dist.tsv", type=Path, required=False,
    help='Path to output file [required]'
)
@click.option(
    '--kmer_size', '-k', default=21, type=str, required=False,
    help='K-mer size to estimate genome distance between all genomes'
)
@click.option(
    '--sketch_size', '-s', default=1000, type=int, required=False,
    help='Sketch size to estimate genome distance between all genomes'
)
def mashdist(fasta, index, output, kmer_size, sketch_size):

    """ Experimental: compute a population graph with NetView and Mash """

    li = LineageIndex(index_file=index)

    with tempfile.TemporaryDirectory() as dirname:
        dirpath = Path(dirname)
        run_cmd(
            f'mash sketch -l -k {kmer_size} -s {sketch_size} '
            f'{fasta} -o {dirpath / "db_lineage"}'
        )

        run_cmd(
            f'mash dist {dirpath / "db_lineage.msh"} {dirpath / "db_lineage.msh"}'
            f' > {dirpath / "db_lineage.tsv"}', shell=True
        )

        df = read_mash_pairwise(dirpath / "db_lineage.tsv")

    matrix = [
        genome_data.dist.values for g2, genome_data in df.groupby('genome2', sort=False)
    ]

    matrix = np.stack(matrix, axis=0)

    np.savetxt(output, matrix, fmt="%.12f", delimiter="\t")


def read_mash_pairwise(file: Path):
    return pandas.read_csv(
        file, sep='\t', usecols=[0, 1, 2], header=None,
        converters={
            0: lambda x: Path(x).stem,
            1: lambda x: Path(x).stem,
        }, names=['genome1', 'genome2', 'dist']
    )
