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
def predict(fastq, sketch, tmp, cores, header, reads):

    """ Online lineage matching and antibiotic susceptibility predictions using
    MinHash from uncorrected nmanopore reads"""

    sketchy = Sketchy()

    try:
        sketchy.predict_nanopore(
            sketch=Path(sketch).resolve(),
            fastq=Path(fastq).resolve(),
            tmp=Path(tmp).resolve(),
            cores=cores,
            header=header,
            nreads=reads,
        )

    except KeyboardInterrupt:
        shutil.rmtree(tmp)
    finally:
        shutil.rmtree(tmp)