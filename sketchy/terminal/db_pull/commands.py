import click
import shutil

from pathlib import Path
from sketchy.cloud.storage import GoogleCloudSketch


@click.command()
@click.option(
    '--path', '--p', default=Path.home() / '.sketchy', type=Path,
    help='Path to sketchy home directory [default: ~/.sketchy ]'
)
def db_pull(path):
    """ Pull sketch databases into local storage. """

    # Managed on GitHub not GCS; keep distinct
    data_path = (path / 'data')
    data_path.mkdir(exist_ok=True, parents=True)

    gcs = GoogleCloudSketch(sketch_path=path / 'db')

    # Data index setup
    gcs.pl.logger.info(f'Copy database index files to {data_path}')
    for file in (Path(__file__).parent.parent.parent / 'data').glob('*.tsv'):
        gcs.pl.logger.debug(f'File check: {file}')
        shutil.copy(file, data_path / file.name)

    # Database sketch setup
    gcs.download()
