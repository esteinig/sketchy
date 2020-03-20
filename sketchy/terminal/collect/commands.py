import click
import pandas

from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns


@click.command()
@click.option(
    '--directory', '-d', type=Path, help='Path to directory to collect {prefix}.data.tsv outputs from'
)
@click.option(
    '--nextflow', '-n', is_flag=True, help='Predictions are from Nextflow pipeline and conform to naming scheme'
)
@click.option(
    '--output', '-o', type=Path, default="summary", help='Output directory for summary files.'
)
@click.option(
    '--prefix', '-p', type=str, default="summary", help='Prefix for summary files.'
)
def collect(directory, nextflow, output, prefix):

    """ Collect predictions and summarize results """

    data = []
    columns = []
    additional = []
    for file in directory.glob("*.data.tsv"):
        name = file.name.replace(".data.tsv", "")

        if nextflow:
            if 'time' in name:
                continue

            prefix_data = name.split(".")

            sketch = directory.name
            reads = prefix_data[-1]
            ranks = prefix_data[-2]
            name = ".".join(
                prefix_data[:-2]
            )

            additional.append(
                pandas.Series((ranks, reads, sketch))
            )

        data.append(
            pandas.read_csv(file, sep='\t', header=0)
        )
        columns.append(name)

    if data:
        all_data = pandas.concat(data, axis=1)
        features = all_data[['feature']]
        new_index = features.iloc[:, 0]

        predictions = all_data[['prediction']]
        predictions.index = new_index
        predictions.columns = columns

        stability = all_data[['stability']]
        stability.index = new_index
        stability.columns = columns

        if additional:
            additional_data = pandas.concat(additional, axis=1)
            additional_data.columns = columns
            additional_data.index = ['ranks', 'reads', 'sketch']

            predictions = pandas.concat([additional_data, predictions])
            stability = pandas.concat([additional_data, stability])

        predictions.T.sort_index().to_csv(f'{prefix}.predictions.tsv', sep='\t')
        stability.T.sort_index().to_csv(f'{prefix}.stability.tsv', sep='\t')
    else:
        print('No files found.')
        exit(1)




