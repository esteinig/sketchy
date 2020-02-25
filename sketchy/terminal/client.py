import click

from .utils import utils
from .list import list
from .pull import pull
from .feature import feature
from .online import online
from .survey import survey
from .run import run
from .plot import plot

VERSION = '0.4.2'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """Sketchy: online lineage and genotype matching using MinHash"""
    pass


terminal_client.add_command(run)
terminal_client.add_command(plot)
terminal_client.add_command(survey)
terminal_client.add_command(pull)
terminal_client.add_command(list)
terminal_client.add_command(online)
terminal_client.add_command(utils)
terminal_client.add_command(feature)
