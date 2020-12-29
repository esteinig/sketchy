import click

from .create import create
from .inspect import inspect

@click.group()
def db():
    """ Reference database construction and inspection """
    pass


online.add_command(monitor)
online.add_command(simulate)
online.add_command(watch)
