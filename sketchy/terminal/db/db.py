import click

from .list import list
from .pull import pull


VERSION = "0.1"


@click.group()
@click.version_option(version=VERSION)
def db():
    """ Pull and list reference data from GCS """
    pass


db.add_command(list)
db.add_command(pull)
