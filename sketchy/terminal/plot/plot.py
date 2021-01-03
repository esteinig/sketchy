import click

from .sssh_diagnostics import sssh_diagnostics
from .genotype_heatmap import evaluation_heatmap


@click.group()
def plot():
    """ Tasks for diagnostics plots and heatmaps """
    pass


plot.add_command(sssh_diagnostics)
plot.add_command(evaluation_heatmap)