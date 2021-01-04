import click

from .sssh_diagnostics import sssh_diagnostics
from .genotype_heatmap import genotype_heatmap
from .ssh_heatmap import ssh_heatmap

@click.group()
def plot():
    """ Tasks for diagnostics plots and heatmaps """
    pass


plot.add_command(ssh_heatmap)
plot.add_command(sssh_diagnostics)
plot.add_command(genotype_heatmap)