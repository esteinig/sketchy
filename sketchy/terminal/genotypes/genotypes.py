import click

from .lineage import lineage
from .prepare import prepare
from .merge import merge
from .drop import drop

VERSION = "0.4.4"


@click.group()
@click.version_option(version=VERSION)
def genotypes():
    """ Tasks for: genotypes index preparation and population mapping """
    pass


genotypes.add_command(merge)
genotypes.add_command(drop)
genotypes.add_command(lineage)
genotypes.add_command(prepare)