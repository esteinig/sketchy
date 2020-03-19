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
def plot_summary(directory):

    """ Generate summary report plots from multiple Sketchy predictions """

    for file in directory.glob("*.data.tsv"):
        prefix = file.name.replace(".data.tsv", "")
        df = pandas.read_csv(file, sep='\t')
