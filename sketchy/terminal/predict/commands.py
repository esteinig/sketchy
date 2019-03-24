import click
import shutil

from pathlib import Path

from sketchy.sketchy import Sketchy


@click.command()
@click.option(
    '--fastq', '-f', help='Input FASTQ file to slice and predict'
)
@click.option(
    '--sketch', '-s', default='saureus-v1.msh', help='MASH sketch to query'
)
@click.option(
    '--data', '-d', help='Index data file to pull genotypes and other data from.'
)
@click.option(
    '--tmp', '-t', default=Path().home() / '.sketchy' / 'tmp',
    help='Temporary dir for slicing Fastq file'
)
@click.option(
    '--cores', '-c', default=8, help='Number of processors for `mash dist`'
)
@click.option(
    '--header', '-h', is_flag=True, help='Print header to online mode STDOUT.'
)
@click.option(
    '--reads', '-r', default=50, help='Number of reads to type.'
)
def predict(fastq, sketch, data, tmp, cores, header, reads):

    """ Online lineage matching from uncorrected nanopore reads"""

    sketchy = Sketchy()

    fastq_path = Path(fastq).resolve()
    sketch_path = Path(sketch).resolve()

    if not fastq_path.exists():
        print(f'File {fastq_path} does not exist.')
        exit(1)

    if not sketch_path.exists():
        print(f'Mash sketch {sketch_path} does not exist.')
        exit(1)

    try:
        sketchy.predict_nanopore(
            sketch=sketch_path,
            fastq=fastq_path,
            tmp=Path(tmp).resolve(),
            cores=cores,
            header=header,
            nreads=reads,
            data=data,
        )

    except KeyboardInterrupt:
        shutil.rmtree(tmp)
    else:
        shutil.rmtree(tmp)
