import click
import pandas
import shutil

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
    '--outdir', '-o', default='sketchy_eval', help='Output directory for evaluation data.'
)
@click.option(
    '--bootstrap', '-b', default=100, help='Number of bootstrap replicates.'
)
@click.option(
    '--reads', '-r', default=100, help='Number of reads to sample for bootstrap replicates and predictions.'
)
@click.option(
    '--shuffle', is_flag=True, help='Shuffle reads before sampling bootstrap replicates.'
)
@click.option(
    '--cores', '-c', default=4, help='Number of processors to run predictions with.'
)
@click.option(
    '--prefix', '-p', default="bootstraps", help='Prefix for output files.'
)
def evaluate(fastq, outdir, bootstrap, reads, sketch, data, shuffle, cores, prefix):
    """ Evaluate a set of FASTQ files with Sketchy """

    outdir = Path(outdir)

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
        bsfiles = evaluator.bootstrap(
            fastq=fq,
            nbootstrap=bootstrap,
            sample_reads=reads,
            shuffle=shuffle
        )
        results = evaluator.predict_bootstraps(
            bsfiles=bsfiles,
            sketch=Path(sketch).resolve(),
            data=Path(data).resolve(),
            cores=cores,
            reads=reads
        )

        results['file'] = [fq.name for _ in results.bootstrap]
        file_data.append(results)

    df = pandas.concat(file_data)

    for d in outdir.glob("*"):
        if d.is_dir():
            shutil.rmtree(d)

    df.to_csv(outdir / f'{prefix}.tab', header=True, index=True, sep='\t')

