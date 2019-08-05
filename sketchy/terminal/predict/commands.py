import click
import shutil
import pysam
import matplotlib.pyplot as plt

from pathlib import Path

from sketchy.evaluation import SampleEvaluator
from sketchy.minhash import MashScore
from sketchy.utils import PoreLogger


@click.command()
@click.option(
    '--fastq', '-f', required=True, type=Path,
    help='Input FASTQ file to predict lineage and traits from.',
)
@click.option(
    '--outdir', '-o', required=True, type=Path,
    help='Output directory for sum of shared hashes data and plots.',
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
    '--cores', '-c', default=2, help='Number of processors for MASH'
)
@click.option(
    '--ncpu', default=4, type=int,
    help='Spin MASH computations into threads, then compute sum'
         ' of shared hashes; not for online compute yet.'
)
@click.option(
    '--top', default=50,  type=int,
    help='Collect the top ranked genome hits by sum of shared hashes to plot.'
)
@click.option(
    '--lineages', default=5,  type=int,
    help='Collect the top ranked lineages aggregated by sum of sums of '
         'shared hashes from the --top collected genomes.'
)
@click.option(
    '--sketchy',  default=Path.home() / '.sketchy', type=Path,
    help='Path to Sketchy home directory [ ~/.sketchy/ ]'
)
@click.option(
    '--data', '-d', type=Path,
    help='Index data file for pull genotypes; optional if template '
         'sketch provided'
)
def predict(
    fastq,
    sketch,
    data,
    outdir,
    tmp,
    keep,
    cores,
    ncpu,
    reads,
    top,
    lineages,
    sketchy
):

    """ Lineage hashing from uncorrected nanopore reads (offline) """

    pl = PoreLogger()
    sketch_path = Path(sketch)
    tmp.mkdir(parents=True, exist_ok=True)

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

    if sketch in ('kleb', 'mrsa', 'tb'):
        sketch_path = sketchy / 'db' / f'{sketch}.default.msh'
        data = sketchy / 'data' / f'{sketch}.data.tsv'

    if not fastq.exists():
        click.echo(f'File {fastq} does not exist.')
        exit(1)

    if not sketch_path.exists():
        click.echo(f'Mash sketch {sketch_path} does not exist.')
        exit(1)

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
            ncpu=ncpu,
        )

        se = SampleEvaluator(
            indir=tmp,
            outdir=outdir,
            limit=reads,
            top=top,
            sketch_data=data
        )

        fig, (ax1, ax2) = plt.subplots(
            nrows=1, ncols=2, figsize=(21.0, 7.0)
        )
        fig.subplots_adjust(hspace=0.5)
        fig.suptitle(f'{tmp.name}')

        se.create_lineage_hitmap(top=lineages, ax=ax1)
        se.create_lineage_plot(top=lineages, ax=ax2)

        plt.tight_layout()

        fig.savefig(
            outdir / 'lineage_plots.pdf',
        )

        se.top_ssh.to_csv(
            outdir / 'lineage_data.tsv', sep='\t'
        )

    except KeyboardInterrupt:
        if not keep:
            shutil.rmtree(tmp)
        exit(0)
    else:
        if not keep:
            shutil.rmtree(tmp)
