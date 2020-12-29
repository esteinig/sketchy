import click

from .utils import utils
from .list import list
from .pull import pull
from .create import create
from .online import online
from .survey import survey
from .run import run
from .plot import plot
from .collect import collect

VERSION = '0.5.0'

@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """Sketchy: online db_lineage and genotype matching using MinHash"""
    pass


terminal_client.add_command(collect)
terminal_client.add_command(run)
terminal_client.add_command(plot)
terminal_client.add_command(survey)
terminal_client.add_command(pull)
terminal_client.add_command(list)
terminal_client.add_command(online)
terminal_client.add_command(utils)
terminal_client.add_command(create)
