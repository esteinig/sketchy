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
    '--outdir', '-o', default='sketchy_eval', help='Output directory for evaluation data.'
)
@click.option(
    '--bootstrap', '-b', default=100, help='Number of bootstrap replicates.'
)
@click.option(
    '--reads', '-r', default=100, help='Number of reads to sample for bootstrap replicates and predictions.'
)
def boot(fastq, outdir, bootstrap, reads):
    """ Create bootstrap replicates for a set of read files (Nextflow)"""

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

    for fq in files:
        _ = evaluator.bootstrap(
            fastq=fq,
            nbootstrap=bootstrap,
            sample_reads=reads,
            shuffle=True
        )

