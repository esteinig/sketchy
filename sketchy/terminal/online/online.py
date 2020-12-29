import click

from .simulate import simulate
from .watch import watch
from .monitor import monitor


@click.group()
def online():
    """ Tasks for: online sequence run monitoring and simulation """
    pass


online.add_command(monitor)
online.add_command(simulate)
online.add_command(watch)
