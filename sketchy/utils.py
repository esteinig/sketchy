import sys
import shlex
import pandas
import pyfastx
import logging
import subprocess

from dateutil import parser
from pathlib import Path
from pysam import FastxFile
from click import Option, UsageError


class PoreLogger:

    def __init__(self, level=logging.ERROR, file: Path = None, name: str = 'Sketchy'):

        logging.basicConfig(
            level=level,
            format="[%(asctime)s] [%(name)-7s]     %(message)s",
            datefmt='%H:%M:%S',
            filename=file
        )

        self.logger = logging.getLogger(name)

class SketchySimulator:

    """ Simulator utility to run a live read simulation / rerun for Sketchy """

    def __init__(
        self,
        fastx: Path,
        fastx_index: Path = None,
        reads_per_file: int = 500,
        barcodes: list = None,
    ):
        self.fastx = fastx
        self.fai, self.build_read = create_fastx_index(fastx)

        # TODO: Remove when comments are fixed in Pyfastx: since v0.6

        if fastx_index:
            self.fastx_index = pandas.read_csv(fastx_index, sep='\t')
        else:
            self.fastx_index = None

        self.reads_per_file = reads_per_file
        self.barcodes = barcodes

    def get_run_index(self, fout: bool = False, sort_by: str = 'start_time'):

        self.fastx_index = pandas.DataFrame(
            [extract_read_data(read) for read in FastxFile(self.fastx)]
        )

        if sort_by:
            self.fastx_index = self.fastx_index.sort_values(sort_by)

        if fout:
            self.fastx_index.to_csv(
                fout, sep='\t', index=False
            )

        return self.fastx_index

    @staticmethod
    def create_header_comment(row):

        """ Recreate the header comment from an index row """

        comment = ''
        for name in row.index.values:
            if name != 'name' and row[name] is not None:
                comment += f'{name}={row[name]} '

        return comment

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

    if callback and output is not None:
        return callback(output)

    return output


def get_output_handle(fpath: str, fastx: bool = False, out: bool = True):

    if fpath == "-":
        if out:
            handle = sys.stdout
        else:
            handle = sys.stdin
    else:
        p = Path(fpath)
        if not p.parent.is_dir():
            raise NotADirectoryError(
                "Directory specified for output file does not exist: {}".format(
                    p.parent
                )
            )

        if fastx:
            handle = FastxFile(p)
        else:
            handle = p.open("w")

    return handle


def create_fastx_index(fastx):

    if is_fasta(fastx):
        return pyfastx.Fasta(
            str(fastx), build_index=True
        ), build_read_fasta
    elif is_fastq(fastx):
        return pyfastx.Fastq(
            str(fastx), build_index=True
        ), build_read_fastq
    else:
        raise ValueError(
            f'Could not determine input file format: {fastx}'
        )


def build_read_fasta(read, comment: str = None):

    """ Build fasta string from pyfastx read """

    return f">{read.name}{' '+comment if comment else ''}\n{read.seq}"


def build_read_fastq(read, comment: str = None):

    """ Build fastq read string from pyfastx read """

    return f"@{read.name}{' '+comment if comment else ''}" \
        f"\n{read.seq}\n+\n{read.qual}"


def is_fasta(fastx: Path):
    with fastx.open() as fin:
        return fin.readline().startswith('>')


def is_fastq(fastx: Path):
    with fastx.open() as fin:
        return fin.readline().startswith('@')


def get_files(path: Path, patterns: list, names: list = None) -> list:

    """ Search a path for a list of patterns and return file paths """

    fnames = []
    for pattern in patterns:
        for file in path.glob(pattern=pattern):
            fname = file.name.replace(
                pattern.replace('*', ''), ''
            )
            if names:
                if fname in names:
                    fnames.append(fname)
                    print(file)
            else:
                fnames.append(fname)
                print(file)

    return fnames


def extract_read_data(read) -> dict:

    """ Extract header fields from pysam.Read """

    data = dict(
        runid=None,
        sampleid=None,
        start_time=None,
        barcode=None,
        ch=None,
        read=None
    )

    if read.comment is None:
        return data

    read_data = read.comment.split(' ')

    for d in read_data:
        for key in data.keys():
            if key in d:
                data[key] = d.replace(f'{key}=', '')

    if data['start_time'] is not None:
        data['start_time'] = parser.parse(
            data['start_time']
        ).strftime('%Y-%m-%dT%H:%M:%SZ')

    data['name'] = read.name

    return data


class MutuallyExclusiveOption(Option):
    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(
            kwargs.pop('mutually_exclusive', [])
        )
        if self.mutually_exclusive:
            kwargs['help'] = kwargs.get('help', '') + (
                ' NOTE: This argument is mutually exclusive with '
                ' arguments: [' + ', '.join(self.mutually_exclusive) + '].'
            )
        super(MutuallyExclusiveOption, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.mutually_exclusive.intersection(opts) and self.name in opts:
            raise UsageError(
                "Illegal usage: `{}` is mutually exclusive with "
                "arguments `{}`.".format(
                    self.name, ', '.join(self.mutually_exclusive)
                )
            )

        return super(MutuallyExclusiveOption, self) \
            .handle_parse_result(ctx, opts, args)
