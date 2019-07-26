import click

from .sort import sort
from .link import link
from .merge import merge
from .watch import watch
from .sketch import sketch
from .survey import survey
from .predict import predict
from .evaluate import evaluate
from .db_list import db_list
from .db_pull import db_pull
from .summary import summary
from .select_fastq import select_fastq

VERSION = '0.3'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """Sketchy: online lineage and genotype matching using MinHash"""
    pass


terminal_client.add_command(sum)
terminal_client.add_command(sort)
terminal_client.add_command(link)
terminal_client.add_command(merge)
terminal_client.add_command(watch)
terminal_client.add_command(survey)
terminal_client.add_command(sketch)
terminal_client.add_command(predict)
terminal_client.add_command(evaluate)
terminal_client.add_command(select_fastq)
terminal_client.add_command(db_list)
terminal_client.add_command(db_pull)