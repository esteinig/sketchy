import click
import shutil

from pathlib import Path

from sketchy.minhash import MashScore


@click.command()
@click.option(
    '--fastq', '-f', required=True, type=Path,
    help='Input FASTQ file to predict from.',
)
@click.option(
    '--sketch', '-s',  type=str, default=None, required=True,
    help='MASH sketch to query or one of the templates: kleb, mrsa, tb'
)
@click.option(
    '--data', '-d', help='Index data file for pull genotypes; '
                         'optional if template sketch', type=Path
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
         ' if true lineage / genotype is known'
)
@click.option(
    '--cores', '-c', default=8, help='Number of processors for MASH'
)
@click.option(
    '--ncpu', default=0, help='Do not compute SSH online; spread over CPUs.'
)
@click.option(
    '--mode', type=str, default="single",
    help='Analysis mode; single, cumulative, direct.'
)
@click.option(
    '--output', '-o', default=Path().cwd() / 'shared_hashes.csv',
    help='Output summary files, use --keep to; deprecated.'
)
@click.option(
    '--show', default=3, help='Show this many lineages in pretty print.'
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
        output,
        show,
        genotype,
        nextflow,
        pretty,
        info
):

    """ Online lineage matching from uncorrected nanopore reads"""

    fastq_path = Path(fastq)
    sketch_path = Path(sketch)

    if sketch in ('kleb', 'mrsa', 'tb'):
        defaults = Path(__file__).parent.parent.parent / 'sketch'
        sketch_path = defaults / f'{sketch}.default.msh'
        data = defaults / f'{sketch}.data.tsv'

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
            top=10,
            mode=mode,
            data=data,
            out=output,
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
