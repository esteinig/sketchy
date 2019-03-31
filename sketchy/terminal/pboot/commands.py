import click
import pandas
import shutil

from sketchy.evaluation import Evaluator
from pathlib import Path


@click.command()

@click.option(
    '--file', '-f', default='sketchy_eval', help='Bootstrap replicate Fastq file.'
)
@click.option(
    '--outdir', '-o', default='sketchy_eval', help='Output directory for evaluation data.'
)
@click.option(
    '--sketch', '-s', default=None, required=True, help='Pathogen sketch to evaluate on.'
)
@click.option(
    '--data', '-d', default=None, required=True, help='Sketch associated data for.'
)
@click.option(
    '--reads', '-r', default=500, required=False, help='Sketch associated data for.'
)
@click.option(
    '--cores', '-c', default=4, help='Number of processors to run predictions with.'
)
@click.option(
    '--prefix', '-p', default="bootstraps", help='Prefix for output files.'
)
def pboot(file, outdir, sketch, reads, data, cores, prefix):
    """ Predict on a bootstrap replicate (Nextflow) """

    outdir = Path(outdir)

    evaluator = Evaluator(
        outdir=outdir
    )

    results = evaluator.predict_bootstraps(
        bsfiles=[Path(file)],
        sketch=Path(sketch).resolve(),
        data=Path(data).resolve(),
        cores=cores,
        reads=reads
    )

    results['bootstrap'] = [Path(file).stem.strip('boot') for _ in results.bootstrap]
    results.to_csv(outdir / f"{prefix}.tab", sep='\t', header=True, index=True)


