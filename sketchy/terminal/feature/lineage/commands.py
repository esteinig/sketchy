import click

from pathlib import Path
from braceexpand import braceexpand
from sketchy.sketchy import LineageIndex
from sketchy.utils import get_files


@click.command()
@click.option(
    '--data', '-d', default=None, type=Path, required=True,
    help='Path to data file to summarize trait data from [required]'
)
@click.option(
    '--lineage', '-l', default="ST93", type=str, required=False,
    help='Trait to show summary for; columns in data file [required]'
)
@click.option(
    '--output', '-o', default=f'lineage.index.tsv', type=Path, required=False,
    help='Path to legacy key file to extract identifiers [lineage.index.tsv]'
)
@click.option(
    '--summary', '-s', is_flag=True,
    help='Print summary of lineage features to terminal [false]'
)
@click.option(
    '--file_path', '-f', default=None, type=Path, required=False,
    help='Path to collect files for this lineage from [none]'
)
@click.option(
    '--pattern', '-p', default='*.fastq.gz', type=str, required=False,
    help='Pattern to match files with their name identifier [*.fastq.gz]'
)
@click.option(
    '--key', '-k', default=None, type=Path, required=False,
    help='Path to legacy key file to extract identifiers [none]'
)
@click.option(
    '--reindex', '-r', is_flag=True,
    help='Reindex the lineage table [none]'
)
def lineage(data, lineage, summary, key, output, file_path, pattern, reindex):

    """ Show a summary of a lineage from the reference sketch data """

    li = LineageIndex(index_file=data)

    if summary:
        df = li.get_summary(lineage)
    else:
        if key:
            df = li.get_key_index(
                lineage=lineage, key_file=key
            )
        else:
            df = li.get_lineage(lineage=lineage)

        if summary:
            li.get_summary(lineage)

        if file_path:
            fnames = get_files(
                path=file_path,
                patterns=braceexpand(pattern),
                names=df.id.tolist() if key else df.uuid.tolist()
            )

            df = df[df.id.isin(fnames) if key else df.uuid.isin(fnames)]

            if key:
                df.id = df.id.astype("category")
                df.id.cat.set_categories(fnames, inplace=True)
                df.sort_values("id", inplace=True)
            else:
                df.uuid = df.uuid.astype("category")
                df.uuid.cat.set_categories(fnames, inplace=True)
                df.sort_values("uuid", inplace=True)

    if reindex:
        df.reset_index(inplace=True)
        if 'idx' in df.columns:
            df.drop(columns='idx', inplace=True)
            df.index.name = 'idx'

    df.to_csv(output, sep='\t', index=True, header=True)