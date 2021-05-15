import click

from .metrics import metrics
from .match import match
from .fastx import fastx
from .database import database
from .online import online
from .survey import survey
from .plot import plot
from .collect import collect

VERSION = '0.5.0'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """Sketchy: genomic neighbor typing using MinHash"""
    pass


terminal_client.add_command(metrics)
terminal_client.add_command(match)
terminal_client.add_command(collect)
terminal_client.add_command(plot)
terminal_client.add_command(survey)
terminal_client.add_command(online)
terminal_client.add_command(fastx)
terminal_client.add_command(database)
