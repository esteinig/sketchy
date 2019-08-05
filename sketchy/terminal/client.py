import click

from .link import link
from .merge import merge
from .watch import watch
from .sketch import sketch
from .survey import survey
from .predict import predict
from .summary import summary
from .concat import concat
from .db_pull import db_pull
from .db_list import db_list
from .select_fastq import select_fastq
from .sort_fastq import sort_fastq

VERSION = '0.3'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """Sketchy: online lineage and genotype matching using MinHash"""
    pass


terminal_client.add_command(summary)
terminal_client.add_command(sort_fastq)
terminal_client.add_command(concat)
terminal_client.add_command(link)
terminal_client.add_command(merge)
terminal_client.add_command(watch)
terminal_client.add_command(survey)
terminal_client.add_command(sketch)
terminal_client.add_command(predict)
terminal_client.add_command(select_fastq)
terminal_client.add_command(db_list)
terminal_client.add_command(db_pull)