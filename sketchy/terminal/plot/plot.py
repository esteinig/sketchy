import click

from .sssh_diagnostics import sssh_diagnostics
from .prediction_heatmap import prediction_heatmap
from .ssh_heatmap import ssh_heatmap

@click.group()
def plot():
    """ Tasks for diagnostics plots and heatmaps """
    pass


plot.add_command(ssh_heatmap)
plot.add_command(sssh_diagnostics)
plot.add_command(prediction_heatmap)