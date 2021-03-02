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
    '--summary', '-s', default=None, type=Path, required=False,
    help='Output file for summary or lineage summary genome table [none]'
)
@click.option(
    '--output', '-o', default=None, type=Path, required=False,
    help='Output lineage summary genome table [none]'
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
    '--id_column', '-id',  default='id', type=str, required=False,
    help='DB identifier column [id]'
)
@click.option(
    '--lineage_column', '-lc',  default='mlst', type=str, required=False,
    help='DB lineage column [mlst]'
)
def inspect(db, lineage, summary, output, file_path, pattern, reindex, id_column, lineage_column):

    """ Interrogate the reference database: summary or subset genotypes """

    li = SketchyDatabase(db_path=db, lineage_column=lineage_column)

    if file_path:
        df = li.get_lineage(lineage=lineage)

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

    if summary:
        df = li.get_summary(lineage)
        df.to_csv(output, sep='\t', index=True, header=True)

    if output:
        df = li.get_lineage(lineage=lineage)
        df.to_csv(output, sep='\t', index=False, header=True)