import click

from .diagnostics import diagnostics
from .evaluation_heatmap import evaluation_heatmap


@click.group()
def plot():
    """ Tasks for diagnostics plots and heatmaps """
    pass


plot.add_command(diagnostics)
plot.add_command(evaluation_heatmap)