import click

from .construct import construct
from .link import link
from .popmap import popmap
from .mashdist import mashdist

VERSION = "0.4.4"


@click.group()
@click.version_option(version=VERSION)
def survey():
    """ Tasks for: creating feature index files from Pathfinder Survey """
    pass


survey.add_command(mashdist)
survey.add_command(popmap)
survey.add_command(construct)
survey.add_command(link)