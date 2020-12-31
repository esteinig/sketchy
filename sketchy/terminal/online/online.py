import click

from .watch import watch
from .monitor import monitor


@click.group()
def online():
    """ Tasks for: online sequence run monitoring and simulation """
    pass


online.add_command(monitor)
online.add_command(watch)
