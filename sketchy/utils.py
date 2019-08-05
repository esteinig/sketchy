import logging
import subprocess
import sys
import shlex
from pathlib import Path
import numpy as np
import scipy.stats


class PoreLogger:

    def __init__(self):

        logging.basicConfig(
            level=logging.INFO,
            format="[%(asctime)s] [%(name)-10s]     %(message)s",
            datefmt='%H:%M:%S',
            filename=None
        )

        self.logger = logging.getLogger(self.__class__.__name__)


BREWER = {
  'blue': [
      "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0",
      "#225ea8", "#253494", "#081d58"
  ],
  'green': [
      "#f7fcb9", "#d9f0a3",	"#addd8e",	"#78c679", "#41ab5d",
      "#238443", "#006837",	"#004529"
  ],
  'red': [
      "#e7e1ef", "#d4b9da", "#c994c7",	"#df65b0", "#e7298a",
      "#ce1256", "#980043", "#67001f"
  ],
  'orange': [
      "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a",
      "#e31a1c", "#bd0026",	"#800026"
  ]
}


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array([v for v in data if v is not None])
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


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
