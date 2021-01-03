import click

from .create import create
from .inspect import inspect


@click.group()
def database():
    """ Reference database construction and inspection """
    pass


database.add_command(create)
database.add_command(inspect)
