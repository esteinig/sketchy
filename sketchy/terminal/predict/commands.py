import click
import shutil

from pathlib import Path

from sketchy.sketchy import Sketchy


@click.command()
@click.option(
    '--fastq', '-f', help='Input FASTQ file to predict from.'
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
@click.option(
    '--single', is_flag=True, help='Single read analysis, instead of cumulative, does not compute score'
)
@click.option(
    '--dist', is_flag=True, help='Use best hits from min hash distance, instead shared hashes.'
)
@click.option(
    '--output', '-o', default=Path().cwd() / 'shared_hashes.csv',
    help='Top shared hash queries per read, not cumulative, does not compute scores.'
)
def predict(fastq, sketch, data, tmp, cores, header, reads, single, output, dist):

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
            fastq=fastq_path,
            sketch=sketch_path,
            data=data,
            tmp=Path(tmp).resolve(),
            cores=cores,
            header=header,
            nreads=reads,
            score=not single,
            out=output,
            sort_by='dist' if dist else 'shared'
        )
    except KeyboardInterrupt:
        shutil.rmtree(tmp)
    else:
        shutil.rmtree(tmp)
