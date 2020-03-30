import click

from .plot_kraken import plot_kraken
from .fx_filter import fx_filter
from .fx_sample import fx_sample
from .fx_sort import fx_sort
from .fx_index import fx_index
from .fx_time import fx_time

VERSION = "0.4.4"


@click.group()
@click.version_option(version=VERSION)
def utils():
    """ Tasks for: supplementary utility functions """
    pass


utils.add_command(fx_filter)
utils.add_command(fx_sample)
utils.add_command(fx_sort)
utils.add_command(fx_index)
utils.add_command(fx_time)
utils.add_command(plot_kraken)