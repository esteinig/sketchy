import click

from .genotype_diagnostics import genotype_diagnostics
from .genotype_heatmap import genotype_heatmap
from .raw_heatmap import raw_heatmap

@click.group()
def plot():
    """ Tasks for diagnostics plots and heatmaps """
    pass


plot.add_command(raw_heatmap)
plot.add_command(genotype_diagnostics)
plot.add_command(genotype_heatmap)