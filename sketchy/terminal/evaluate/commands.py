import click
import pandas

from sketchy.evaluation import Evaluator
from pathlib import Path


@click.command()
@click.option(
    '--fastq', '-f', default=None, required=True, help='Input path glob (e.g. *.fq) for Fastq, or single Fastq.'
)
@click.option(
    '--sketch', '-s', default=None, required=True, help='Pathogen sketch to evaluate on.'
)
@click.option(
    '--data', '-d', default=None, required=True, help='Sketch associated data for.'
)
@click.option(
    '--outdir', '-o', default='sketchy_eval', help='Output directory for evaluations data.'
)
@click.option(
    '--bootstrap', '-b', default=100, help='Number of bootstrap replicates.'
)
@click.option(
    '--reads', '-r', default=1000, help='Number of reads to sample for bootstrap replicates and predictions.'
)
@click.option(
    '--cores', '-c', default=4, help='Number of processors to run predictions with.'
)
def evaluate(fastq, outdir, bootstrap, reads, sketch, data, cores):
    """ Evaluate a set of FASTQ files with Sketchy """

    if '*' in fastq:
        files = list(
            Path().glob(fastq)
        )
    else:
        files = [Path(fastq)]

    evaluator = Evaluator(
        outdir=outdir
    )

    file_data = []
    for fq in files:
        bsfiles = evaluator.bootstrap(fastq=fq, nb=bootstrap, reads=reads)
        results = evaluator.predict_bootstraps(bsfiles=bsfiles, sketch=sketch, data=data, cores=cores)
        results.id = [fq.stem for _ in results.bootstrap]
        file_data.append(results)

    df = pandas.concat(file_data)

    print(df)