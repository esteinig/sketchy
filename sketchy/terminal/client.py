import click

from .sort import sort
from .predict import predict

VERSION = '0.1'


@click.group()
@click.version_option(version=VERSION)
def terminal_client():
    """ Sketchy: online lineage matchting and
    antimicrobial susceptibility predictions using MinHash"""
    pass


terminal_client.add_command(sort)
terminal_client.add_command(predict)