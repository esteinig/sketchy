import click
from pathlib import Path
from sketchy.storage import GoogleCloudSketch


@click.command()
@click.option(
    '--path', '--p', default=Path.home() / '.sketchy', type=Path,
    help='Path to sketchy home directory [~/.sketchy ]'
)
def pull(path):

    """ Pull default sketch collection into local storage """

    path.mkdir(exist_ok=True, parents=True)
    GoogleCloudSketch(sketch_path=path, full=False).pull()
