import click

from pathlib import Path
from braceexpand import braceexpand
from sketchy.sketchy import SketchyDatabase
from sketchy.utils import get_files


@click.command()
@click.option(
    '--db', '-d', default=None, type=Path, required=True,
    help='Path to Sketchy database [required]'
)
@click.option(
    '--lineage', '-l', default="ST93", type=str, required=False,
    help='Lineage to show summary for [required]'
)
@click.option(
    '--output', '-o', default=f'db_lineage.index.tsv', type=Path, required=False,
    help='Path to legacy key file to extract identifiers [lineage.index.tsv]'
)
@click.option(
    '--summary', '-s', is_flag=True,
    help='Print summary of db_lineage features to terminal [false]'
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
    '--reindex', '-r', is_flag=True,
    help='Reindex the lineage table [none]'
)
@click.option(
    '--id_column', '-ic',  default='id', type=str, required=False,
    help='DB identifier column [id]'
)
@click.option(
    '--lineage_column', '-lc',  default='mlst', type=str, required=False,
    help='DB lineage column [mlst]'
)
def inspect(data, lineage, summary, output, file_path, pattern, reindex, id_column, lineage_column):

    """ Interrogate the reference database """

    li = SketchyDatabase(genotype_file=data, lineage_column=lineage_column)

    if summary:
        df = li.get_summary(lineage)
    else:
        df = li.get_lineage(lineage=lineage)

        if summary:
            li.get_summary(lineage)

        if file_path:
            fnames = get_files(
                path=file_path,
                patterns=braceexpand(pattern),
                names=df[id_column].tolist()
            )

            df = df[df[id_column].isin(fnames)]

            df[id_column] = df[id_column].astype("category")
            df[id_column].cat.set_categories(fnames, inplace=True)
            df.sort_values(id_column, inplace=True)

    if reindex:
        df.reset_index(inplace=True)
        if 'idx' in df.columns:
            df.drop(columns='idx', inplace=True)
            df.index.name = 'idx'

    df.to_csv(output, sep='\t', index=True, header=True)