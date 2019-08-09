import click
import shutil
import pysam

from pathlib import Path

from sketchy.evaluation import SampleEvaluator
from sketchy.minhash import MashScore
from sketchy.utils import PoreLogger, get_total_reads


@click.command()
@click.option(
    '--prefix', '-p', '-o', required=True, type=str, default='sketchy',
    help='Output prefix for sum of shared hashes data.',
)
@click.option(
    '--fastq', '-f', required=True, type=Path,
    help='Input FASTQ file to predict lineage and traits from.',
)
@click.option(
    '--sketch', '-s',  type=str, default=None, required=True,
    help='MASH sketch file to query; or a template, one of: kleb, mrsa, tb'
)
@click.option(
    '--reads', '-r', default=1000, help='Number of reads to type.', type=int
)
@click.option(
    '--tmp', '-t', default=Path().cwd() / 'tmp', type=Path,
    help='Temporary directory for sum of shared hashes per read output.'
)
@click.option(
    '--keep', '-k', is_flag=True,
    help='Keep temporary folder with per read shared hashes.'
)
@click.option(
    '--threads', '-t', default=4, type=int,
    help='Spin MASH into threads, compute post-hoc sum '
         'of shared hashes.'
)
@click.option(
    '--top', default=10,  type=int,
    help='Collect the top ranked genome hits by sum of shared hashes to plot.'
)
@click.option(
    '--cores', default=2, help='Number of processors for MASH'
)
@click.option(
    '--data', type=Path,
    help='Index data file for pull genotypes; optional if template '
         'sketch provided'
)
@click.option(
    '--sketchy',  default=Path.home() / '.sketchy', type=Path,
    help='Path to Sketchy home directory [ ~/.sketchy/ ]'
)
def predict(
    fastq,
    sketch,
    data,
    prefix,
    tmp,
    keep,
    cores,
    threads,
    reads,
    top,
    sketchy
):

    """ Lineage hashing from uncorrected nanopore reads (offline) """

    if reads == 0:
        reads = None

    pl = PoreLogger()
    sketch_path = Path(sketch)
    tmp.mkdir(parents=True, exist_ok=True)

    if sketch in ('kleb', 'mrsa', 'tb'):
        sketch_path = sketchy / 'db' / f'{sketch}.default.msh'
        data = sketchy / 'data' / f'{sketch}.data.tsv'

    if not fastq.exists():
        click.echo(f'File {fastq} does not exist.')
        exit(1)

    if not sketch_path.exists():
        click.echo(f'Mash sketch {sketch_path} does not exist.')
        exit(1)

    if fastq.suffix == '.gz':
        # Unpack into temporary directory
        tmp_path = tmp / fastq.with_suffix('')
        pl.logger.debug(f'Decompressing file {fastq} to {tmp_path}')
        with pysam.FastxFile(fastq) as fin, \
                open(tmp / fastq.with_suffix(''), mode='w') as fout:
            for entry in fin:
                string_out = str(entry)
                if not string_out.endswith('\n'):
                    string_out += '\n'
                fout.write(string_out)

        fastq = tmp_path

    if reads is None:
        reads = get_total_reads(fastq)

    try:

        ms = MashScore()
        pl.logger.info('Compute min-wise shared hashes against sketch ...')

        _ = ms.run(
            fastq=fastq,
            nreads=reads,
            sketch=sketch_path,
            cores=cores,
            top=top,  # direct mode only
            mode='single',
            data=data,
            tmpdir=tmp,
            ncpu=threads,
        )

        se = SampleEvaluator(
            indir=tmp,
            limit=reads,
            top=top,
            sketch_data=data
        )

        se.top_ssh.to_csv(
            f'{str(prefix) + ".tsv"}', sep='\t'
        )

    except KeyboardInterrupt:
        if not keep:
            shutil.rmtree(tmp)
        exit(0)
    else:
        if not keep:
            shutil.rmtree(tmp)
