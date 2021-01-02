import click

from .plot_stream import plot_stream
from .plot_heatmap import plot_heatmap

@click.group()
def diagnostics():
    """ Tasks for diagnostics and plots """
    pass


diagnostics.add_command(plot_stream)
diagnostics.add_command(plot_heatmap)