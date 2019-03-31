import click


from .sort import sort
from .boot import boot
from .psum import psum
from .pplot import pplot
from .pboot import pboot
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
terminal_client.add_command(boot)
terminal_client.add_command(psum)
terminal_client.add_command(pplot)
terminal_client.add_command(pboot)
terminal_client.add_command(sketch)
terminal_client.add_command(predict)
terminal_client.add_command(evaluate)
