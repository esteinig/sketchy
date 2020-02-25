import click

from .lineage import lineage
from .prepare import prepare
from .merge import merge
from .drop import drop

VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def feature():
    """ Tasks for: feature index preparation and population mapping """
    pass


feature.add_command(merge)
feature.add_command(drop)
feature.add_command(lineage)
feature.add_command(prepare)