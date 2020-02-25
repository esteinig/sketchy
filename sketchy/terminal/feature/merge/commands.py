import click
import os
import pandas
import logging

from pathlib import Path
from sketchy.utils import run_cmd, PoreLogger


@click.command()
@click.option(
    '--sketch', '-s', type=Path, required=True,
    help='Path to sketch file to parse indices from [required]'
)
@click.option(
    '--features', '-f', type=Path, required=True,
    help='Path to genotype feature file to merge indices with sketch [required]'
)
@click.option(
    '--key', '-k', type=Path, required=False,
    help='Legacy key file to translate UUIDs [dep.]'
)
@click.option(
    '--index_column', '-i', type=str, default='idx',
    help='Feature index column to merge indices on [idx]'
)
@click.option(
    '--prefix', '-p', type=str, default='sketchy.info',
    help='Prefix for output file: {prefix}.tsv [sketchy]'
)
@click.option(
    '--verbose', '-v', is_flag=True,
    help='Enable verbose output for merge operations'
)
def merge(sketch, features, key, prefix, index_column, verbose):

    """ Merge sketch and feature data by common indices """

    pl = PoreLogger(level=logging.INFO if verbose else logging.ERROR).logger

    pl.info(f'Extracting data from sketch: {sketch}')
    run_cmd(f'mash info -t {sketch} > {prefix}.mashinfo', shell=True)

    pl.info(f'Reading and converting data indices from sketch')
    converters = {'id': lambda x: Path(x).stem}
    mash_info = pandas.read_csv(
        f'{prefix}.mashinfo',
        sep='\t',
        header=None,
        skiprows=1,
        index_col=0,
        engine='c',
        usecols=[2],
        names=['id'],
        converters=converters,
    )

    pl.info(f'Assigning sequential indices to index column: `idx`')
    mash_info['idx'] = [i for i in range(len(mash_info))]
    mash_info['ids'] = mash_info.index.tolist()

    nsketch = len(mash_info)

    pl.info(f'Ordered merge on column {index_column} with feature file {features}')
    d = pandas.read_csv(features, sep='\t')

    ndata = len(d)

    print(mash_info)
    print(d)

    mash_info = d.merge(
        mash_info, left_on=index_column, right_on='idx', how='inner'
    )
    pl.info('Merged data and sketch information')

    if key is not None:
        key_table = pandas.read_csv(
            key, sep='\t', header=0
        )
        mash_info = mash_info.merge(
            key_table, left_on='ids', right_on='uuid'
        )
        mash_info.drop(columns=['uuid', 'fasta'], inplace=True)
        mash_info.rename(columns={'id': 'key'}, inplace=True)

    print(mash_info)
    pl.info(f'Writing merged feature index to: {prefix}.tsv')
    mash_info.to_csv(
        f'{prefix}.tsv',
        sep='\t',
        header=True,
        index=False,
    )

    pl.info(f'Merged sketch data ({nsketch}) and feature data ({ndata})')
    pl.info(f'Final sketch and feature size is {len(mash_info)}')
    pl.info(f'Removed features not present in sketch: {len(mash_info) - ndata}')
    pl.info(f'Removing temporary file {prefix}.mashinfo')
    os.remove(f'{prefix}.mashinfo')
