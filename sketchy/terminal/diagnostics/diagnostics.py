import click

from .plot_stream import plot_stream


@click.group()
def diagnostics():
    """ Tasks for diagnostics and plots """
    pass


diagnostics.add_command(plot_stream)
