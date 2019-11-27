import click
import pandas
import pysam
import uuid

from pathlib import Path

from sketchy.sketchy import Sketchy


@click.command()
@click.option(
    '--fastx', '-f', required=True, type=Path,
    help='Fast{a,q} file for prediction with Sketchy',
)
@click.option(
    '--sketch', '-s', required=True, type=Path,
    help='Reference sketch for prediction with Sketchy',
)
@click.option(
    '--reads', '-r', default=None, type=int,
    help='Number of reads to predict from with Sketchy [1000]'
)
@click.option(
    '--threads', '-t', default=4, type=int,
    help='Number of threads for processing with MASH [4]'
)
@click.option(
    '--prefix', '-p', default='sketchy', type=str,
    help='Prefix for ranked sum of shared hashes output with Sketchy [sketchy]'
)
@click.option(
    '--memory', is_flag=True,
    help='Memory efficient computation; significantly slower [false]'
)
@click.option(
    '--ranks', default=10, type=int,
    help='Top ranking sum of shared hashes to extract for evaluation [10]'
)
def test(
    fastx,
    reads,
    sketch,
    ranks,
    prefix,
    memory,
    threads,
):

    """ Task for running prototype code for Sketchy """

    tmp = Path(prefix + '_' + str(
        uuid.uuid4()
    ))

    sketchy = Sketchy(
        fastx=fastx, sketch=sketch, reads=reads, tmp=tmp, memory=memory
    )

    read_list = sketchy.extract_reads()

    if not (tmp / 'mash.tsv').exists():
        results = sketchy.compute_shared_hashes(
            query=read_list, threads=threads,
        )
    else:
        results = tmp / 'mash.tsv'

    ssh = sketchy.compute_ssh(file=results, ranks=ranks)

    ssh.to_csv(f'{prefix}.tsv', sep='\t')
