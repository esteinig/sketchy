import click
import pandas
import pysam

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
    '--reads', '-r', default=1000, type=int,
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
    '--ranks', default=10, type=int,
    help='Top ranking sum of shared hashes to extract for evaluation [10]'
)
@click.option(
    '--p_value', default=1e-06, type=float,
    help='P-value to pre-filter output of shared hashes from MASH [1e-06]'
)
@click.option(
    '--tmp', default=Path().cwd() / 'tmp', type=Path,
    help='Temporary directory for intermediary processing files [tmp]'
)
def test(
    fastx,
    reads,
    sketch,
    ranks,
    prefix,
    p_value,
    threads,
    tmp
):

    """ Task for running prototype code for Sketchy """

    sketchy = Sketchy(
        fastx=fastx,  reads=reads, tmp=tmp
    )

    read_list = sketchy.extract_reads()

    results = sketchy.compute_shared_hashes(
        reference=sketch,
        query=read_list,
        p_value=p_value,
        threads=threads,
    )

    ssh = sketchy.compute_ssh(file=results, ranks=ranks)

    ssh.to_csv(f'{prefix}.tsv', sep='\t')
