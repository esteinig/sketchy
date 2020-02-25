import click
import time

from pathlib import Path
from sketchy.watcher import Watcher


@click.command()
@click.option(
    '--directory', '-d', default=None, type=Path, required=True,
    help='Path to directory to watch [required]'
)
@click.option(
    '--regex', '-r', default=r".*\.fastq$", type=str, required=False,
    help='Regex to identify read files [.*\.fastq$]'
)
@click.option(
    '--now', '-n', is_flag=True, help='Disable waiting for file completion.'
)
def watch(directory, regex, now):

    """ Watch a directory and output file paths to pipe into Sketchy """

    regexes = [regex, ]

    watcher = Watcher()
    try:
        watcher.watch_path(
            path=directory, regexes=regexes, wait=not now,
            callback=lambda x: print(x.absolute())
        )
        while True:
            time.sleep(10)
    except KeyboardInterrupt:
        watcher.watcher.stop()
        exit(0)
    except (OSError, SystemError, SystemExit):
        watcher.watcher.stop()
        exit(1)
