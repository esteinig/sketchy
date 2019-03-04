import click
import Path

from sketchy.sketchy import Sketchy

@click.group()
@click.option(
    '--fastq', '-f', help='Input FASTQ file to slice and predict'
)
@click.option(
    '--sketch', '-s', default='saureus-v1.msh', help='MASH sketch to query'
)
@click.option(
    '--tmp', '-t', default=Path().home() / '.sketchy' / 'tmp',
    help='Temporary dir for slicing FASTQ'
)
@click.option(
    '--cores', '-c', default=8, help='Number of processors for `mash dist`'
)
@click.option(
    '--extension', '-e', default='.fq', help='Extension of read files to glob.'
)
@click.option(
    '--header', '-head', help='Print header to online mode STDOUT.'
)
def predict(fastq, sketch, tmp, cores, extension, header):

    sketchy = Sketchy()

    sketchy.predict_nanopore(
        sketch='saureus-v1.msh',
        fastq='st243_r9.4_reads_sorted.fq',
        tmp='./tmp',
        cores=16,
        extension='.fq',
        header=True
    )
