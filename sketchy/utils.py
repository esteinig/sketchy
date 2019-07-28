import logging
import subprocess
import sys
import shlex
from pathlib import Path


class PoreLogger:

    def __init__(self):

        logging.basicConfig(
            level=logging.INFO,
            format="[%(asctime)s] [%(name)-10s]     %(message)s",
            datefmt='%H:%M:%S',
            filename=None
        )

        self.logger = logging.getLogger(self.__class__.__name__)


def run_cmd(cmd, callback=None, watch=False, background=False, shell=False):

    """Runs the given command and gathers the output.

    If a callback is provided, then the output is sent to it, otherwise it
    is just returned.

    Optionally, the output of the command can be "watched" and whenever new
    output is detected, it will be sent to the given `callback`.

    Returns:
        A string containing the output of the command, or None if a `callback`
        was given.
    Raises:
        RuntimeError: When `watch` is True, but no callback is given.

    """
    if watch and not callback:
        raise RuntimeError(
            "You must provide a callback when watching a process."
        )

    output = None
    try:
        if shell:
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        else:
            proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)

        if background:
            # Let task run in background and return pmid for monitoring:
            return proc.pid, proc

        if watch:
            while proc.poll() is None:
                line = proc.stdout.readline()
                if line != "":
                    callback(line)

            # Sometimes the process exits before we have all of the output, so
            # we need to gather the remainder of the output.
            remainder = proc.communicate()[0]
            if remainder:
                callback(remainder)
        else:
            output = proc.communicate()[0]
    except:
        err = str(sys.exc_info()[1]) + "\n"
        output = err

    if callback and output is not None:
        return callback(output)

    return output


def get_default_sketch(index) -> (Path, Path):
    """ Return the default database and index for Sketchy """
    return (Path(__file__).parent / 'sketch' / f'{index}.default.msh',
        Path(__file__).parent / 'sketch' / f'{index}.data.tsv')
