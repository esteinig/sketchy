import click

from .get import get
from .create import create
from .inspect import inspect


@click.group()
def database():
    """ Reference database construction and inspection """
    pass


database.add_command(get)
database.add_command(create)
database.add_command(inspect)
