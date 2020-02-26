import click
from pathlib import Path
from sketchy.storage import GoogleCloudSketch


@click.command()
@click.option(
    '--path', '--p', default=Path.home() / '.sketchy', type=Path,
    help='Path to sketchy home directory [~/.sketchy ]'
)
@click.option(
    '--full', '-f', is_flag=True,
    help='Pull the full default sketch collections [false]'
)
def pull(path, full):

    """ Pull default sketch collection into local storage """

    path.mkdir(exist_ok=True, parents=True)
    GoogleCloudSketch(sketch_path=path, full=full).pull()
