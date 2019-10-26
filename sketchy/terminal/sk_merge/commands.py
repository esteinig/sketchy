import click
import pandas

from pathlib import Path


@click.command()
@click.option(
    '--survey', '-s', help='Input survey file from: sketchy survey', type=Path
)
@click.option(
    '--key', '-k', help='Input key file with UUIDs from: sketchy link', type=Path
)
@click.option(
    '--out', '-o', default=Path('sketchy.merge.tsv'),
    help='Output data index file for Sketchy', type=Path
)
def sk_merge(survey, key, out):
    """ Create Sketchy data by merging survey with a key file  """

    s = pandas.read_csv(survey, '\t', index_col=0)
    s.index.name = 'id'

    k = pandas.read_csv(key, sep='\t', index_col=0)
    k.reset_index(level=0, inplace=True)

    df = s.merge(k, on='id')

    to_drop = ['id']

    if 'fasta' in df.columns:
        to_drop.append('fasta')

    df = df.drop(columns=to_drop)

    df.set_index('uuid').to_csv(out, sep='\t')