import click

from .sort import sort
from .sketch import sketch
from .predict import predict
from .evaluate import evaluate

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """Sketchy: online lineage and genotype matching using MinHash"""
    pass


terminal_client.add_command(sort)
terminal_client.add_command(sketch)
terminal_client.add_command(predict)
terminal_client.add_command(evaluate)
