import click
import logging
import random

from collections import Counter
from pysam import FastxFile
from pathlib import Path
from sketchy.utils import get_output_handle, is_fasta, is_fastq, create_fastx_index


@click.command()
@click.option(
    "--fastx",
    "-f",
    type=Path,
    help="Path to Fast{a,q} input file.",
    required=True,
)
@click.option(
    "--output",
    "-o",
    type=str,
    help="Output to Fast{a,q} file. Default stdout [-]",
    default="-",
)
@click.option(
    "--sample",
    "-s",
    type=int,
    help="Sample size in number of reads [1000].",
    default=None,
)
@click.option(
    "--replacement",
    "-r",
    help="Sample with replacement [false].",
    is_flag=True
)
@click.option(
    "--seed",
    help="Seed for sampling function [system time].",
    type=int,
    default=None
)
def sample(fastx, output, sample, replacement, seed):

    """ Sample reads with or without replacement (samples held in memory) """

    random.seed(seed)

    fx, build_read = create_fastx_index(fastx)

    read_names = [read.name for read in fx]

    if sample is None:
        sample = len(read_names)

    if len(read_names) < sample:
        print(
            f'Error: Sample ({sample}) must be <='
            f' number of reads ({len(read_names)})'
        )
        exit(1)

    if replacement:
        sampled_reads = random.choices(read_names, k=sample)
    else:
        sampled_reads = random.sample(read_names, k=sample)

    with get_output_handle(output) as fout:
        for read_name in sampled_reads:
            read_str = build_read(
                fx[read_name]
            )
            fout.write(read_str + '\n')

