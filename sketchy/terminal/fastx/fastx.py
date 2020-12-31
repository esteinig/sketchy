import click


from .filter import filter
from .sample import sample
from .sort import sort
from .index import index
from .time import time


@click.group()
def fastx():
    """ Tasks for sequence read file manipulation """
    pass


fastx.add_command(filter)
fastx.add_command(sample)
fastx.add_command(sort)
fastx.add_command(index)
fastx.add_command(time)