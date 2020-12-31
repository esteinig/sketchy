import click

from .construct import construct
from .link import link
from .popmap import popmap
from .mashdist import mashdist


@click.group()
def survey():
    """ Tasks for public database surveillance and sketch updates """
    pass


survey.add_command(mashdist)
survey.add_command(popmap)
survey.add_command(construct)
survey.add_command(link)