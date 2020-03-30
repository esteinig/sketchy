import click

from .simulate import simulate
from .watch import watch
from .monitor import monitor

VERSION = "0.4.4"


@click.group()
@click.version_option(version=VERSION)
def online():
    """ Tasks for: online sequence run monitoring and simulation """
    pass


online.add_command(monitor)
online.add_command(simulate)
online.add_command(watch)
