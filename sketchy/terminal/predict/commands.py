import click
import shutil

from pathlib import Path

from sketchy.minhash import MashScore


@click.command()
@click.option(
    '--fastq', '-f', required=True, type=Path,
    help='Input FASTQ file to predict lineage and traits from.',
)
@click.option(
    '--sketch', '-s',  type=str, default=None, required=True,
    help='MASH sketch file to query; or a template, one of: kleb, mrsa, tb'
)
@click.option(
    '--data', '-d', type=Path,
    help='Index data file for pull genotypes; '
         'optional if template sketch provided'
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
    help='Keep temporary folder for: sketchy evaluate'
)
@click.option(
    '--cores', '-c', default=2, help='Number of processors for MASH'
)
@click.option(
    '--ncpu', default=0, type=int,
    help='Do not compute SSH online; spread over CPUs; not online yet.'
)
@click.option(
    '--mode', type=str, default="single",
    help='Analysis mode; single, cumulative, direct.'
)
@click.option(
    '--show', default=10,
    help='Show this many lineages in pretty print; set to 3 for TB.'
)
@click.option(
    '--genotype',  '-g', is_flag=True, help='Show genotype in pretty print.'
)
@click.option(
    '--nextflow',  '-n', is_flag=True,
    help='Disable sequential online computation for '
         'distributed compute with Nextflow.'
)
@click.option(
    '--pretty',  '-p', is_flag=True, help='Pretty print output.'
)
@click.option(
    '--info',  '-i', is_flag=True, help='Read length and timestamp (adds IO).'
)
@click.option(
    '--sketchy',  default=Path.home() / '.sketchy', type=Path,
    help='Path to Sketchy home directory [ ~/.sketchy/ ]'
)
def predict(
        fastq,
        sketch,
        data,
        tmp,
        keep,
        cores,
        ncpu,
        reads,
        mode,
        show,
        genotype,
        nextflow,
        pretty,
        info,
        sketchy
):

    """ Online lineage matching from uncorrected nanopore reads"""

    fastq_path = Path(fastq)
    sketch_path = Path(sketch)

    if sketch in ('kleb', 'mrsa', 'tb'):
        sketch_path = sketchy / 'db' / f'{sketch}.default.msh'
        data = sketchy / 'data' / f'{sketch}.data.tsv'

    if not fastq_path.exists():
        click.echo(f'File {fastq_path} does not exist.')
        exit(1)

    if not sketch_path.exists():
        click.echo(f'Mash sketch {sketch_path} does not exist.')
        exit(1)

    tmp.mkdir(parents=True, exist_ok=True)

    try:

        ms = MashScore()
        _ = ms.run(
            fastq=fastq_path,
            nreads=reads,
            sketch=sketch_path,
            cores=cores,
            top=10,  # direct mode
            mode=mode,
            data=data,
            out=None,
            sort_by='shared',
            tmpdir=tmp,
            show_top=show,
            show_genotype=genotype,
            nextflow=nextflow,
            ncpu=ncpu,
            pretty=pretty,
            info=info
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
