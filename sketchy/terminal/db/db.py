import click

from .create import create
from .inspect import inspect

@click.group()
def db():
    """ Reference database construction and inspection """
    pass


db.add_command(create)
db.add_command(inspect)
