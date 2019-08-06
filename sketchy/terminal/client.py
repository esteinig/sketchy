import click

from .plot import plot
from .sk_link import sk_link
from .sk_merge import sk_merge
from .sk_sketch import sk_sketch
from .sk_survey import sk_survey
from .predict import predict
from .sk_cat import sk_cat
from .db_pull import db_pull
from .db_list import db_list
from .fq_filter import fq_filter
from .fq_sample import fq_sample
from .fq_sort import fq_sort

VERSION = '0.3'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """Sketchy: online lineage and genotype matching using MinHash"""
    pass


terminal_client.add_command(plot)
terminal_client.add_command(fq_sort)
terminal_client.add_command(fq_sample)
terminal_client.add_command(sk_cat)
terminal_client.add_command(sk_link)
terminal_client.add_command(sk_merge)
terminal_client.add_command(sk_survey)
terminal_client.add_command(sk_sketch)
terminal_client.add_command(predict)
terminal_client.add_command(fq_filter)
terminal_client.add_command(db_list)
terminal_client.add_command(db_pull)