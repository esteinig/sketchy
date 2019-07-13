import click
import shutil

from pathlib import Path

from sketchy.minhash import MashScore


@click.command()
@click.option(
    '--fastq', '-f', help='Input FASTQ file to predict from.', type=Path,
)
@click.option(
    '--sketch', '-s', default='saureus-v1.msh', help='MASH sketch to query', type=Path
)
@click.option(
    '--data', '-d', help='Index data file for pull genotypes.', type=Path
)
@click.option(
    '--reads', '-r', default=500, help='Number of reads to type.', type=int
)
@click.option(
    '--tmp', '-t', default=Path().cwd() / 'tmp', type=Path,
    help='Temporary directory for slicing read file.'
)
@click.option(
    '--keep', '-k', is_flag=True,
    help='Keep temporary folder with match tables for each read / read set.'
)
@click.option(
    '--cores', '-c', default=8, help='Number of processors for `mash dist`'
)
@click.option(
    '--mode', type=str, default="single",
    help='Analysis mode; single, cumulative, direct.'
)
@click.option(
    '--dist', is_flag=True,
    help='Use smallest MinHash distance, instead of most shared hashes.'
)
@click.option(
    '--output', '-o', default=Path().cwd() / 'shared_hashes.csv',
    help='Output summary files, use --keep to'
)
def predict(
        fastq,
        sketch,
        data,
        tmp,
        keep,
        cores,
        reads,
        mode,
        output,
        dist
):

    """ Online lineage matching from uncorrected nanopore reads"""

    fastq_path = Path(fastq)
    sketch_path = Path(sketch)

    if not fastq_path.exists():
        print(f'File {fastq_path} does not exist.')
        exit(1)

    if not sketch_path.exists():
        print(f'Mash sketch {sketch_path} does not exist.')
        exit(1)

    tmp.mkdir(parents=True, exist_ok=True)

    try:

        ms = MashScore()
        _ = ms.run(
            fastq=fastq_path,
            nreads=reads,
            sketch=sketch_path,
            cores=cores,
            top=10,
            mode=mode,
            data=data,
            out=output,
            sort_by='dist' if dist else 'shared',
            tmpdir=tmp,
            online=False
        )

    except KeyboardInterrupt:
        if not keep:
            shutil.rmtree(tmp)
    except AttributeError:
        # KeyboardInterrupt SIGKILL to running Popen.PIPE
        if not keep:
            shutil.rmtree(tmp)
    else:
        if not keep:
            shutil.rmtree(tmp)
