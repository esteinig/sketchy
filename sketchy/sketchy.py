"""
================================
Access point module for Sketchy
================================

TODO: do not jointly count genotypes, count each individually!

"""

import re
import pandas
import random
import dateutil

from tqdm import tqdm

from pathlib import Path

from sketchy.utils import run_cmd

from Bio import SeqIO
from colorama import Fore

Y = Fore.YELLOW
R = Fore.RED
G = Fore.GREEN
C = Fore.CYAN

LR = Fore.LIGHTRED_EX
LC = Fore.LIGHTCYAN_EX
LB = Fore.LIGHTBLUE_EX
LG = Fore.LIGHTGREEN_EX
LY = Fore.LIGHTYELLOW_EX
LM = Fore.LIGHTMAGENTA_EX

RE = Fore.RESET


class Sketchy:
    """ Main access interface to Sketchy """

    def __init__(
            self,
            survey_result: Path or str = None
    ):

        pass

    def create_refdb_sketch(
            self,
            domains: list = None
        ):

        pass

    @staticmethod
    def sort_fastq(
        file: str = None,
        fastq: str = None,
        shuffle: bool = False,
        regex='start_time=(.*)Z',
        first_read: bool = False,
        outfile: Path = None,
    ) -> list or str:
        """ Not very efficient FASTQ sort and shuffle """

        dates = []
        ids = []
        lengths = []
        records = {}

        out = run_cmd(f'wc -l {file}').decode("utf-8").strip().split()

        try:
            total = int(out[0]) // 4
        except (IndexError, TypeError):
            total = None

        with open(file, "r") as input_handle:
            for record in tqdm(
                    SeqIO.parse(input_handle, "fastq"), total=total
            ):
                try:
                    # Extract start time
                    timestr = re.search(regex, record.description)
                    time = timestr.group(1).strip().replace('start_time=', '')
                    dtime = dateutil.parser.parse(time)
                except AttributeError:
                    dtime = None

                dates.append(dtime)
                ids.append(record.id)
                lengths.append(
                    len(record.seq)
                )

                if not first_read:
                    # Do not store the records if
                    # detect first read time is on
                    records[record.id] = record
                #
                # if not fastq:
                #     print(
                #         record.id, time, len(record.seq)
                #     )

        df = pandas.DataFrame(
            data={
                'date': dates,
                'read': ids,
                'length': lengths
            }
        ).sort_values(by='date')

        pandas.set_option('display.max_rows', 120)

        if first_read:
            print(
                df
            )
            if outfile:
                df.to_csv(outfile, sep='\t')
        else:
            df = df.reset_index()

            recs = [
                records[read] for read in df['read']
            ]

            if shuffle:
                random.shuffle(recs)

            with open(fastq, "w") as output_handle:
                SeqIO.write(recs, output_handle, 'fastq')

        return fastq

    @staticmethod
    def _slice_fastq(
        fastq: Path,
        nreads: int = 1,
        tmpdir: Path = Path().cwd() / 'tmp',
        cumulative: bool = True,
    ) -> Path:

        if cumulative:
            fpath = tmpdir / f'reads_{nreads}.fq'
            nreads = 4*nreads
            run_cmd(
                f'head -n {nreads} {fastq} > {fpath}', shell=True
            )
        else:
            fpath = []
            with fastq.open('r') as input_handle:
                for record in SeqIO.parse(input_handle, "fastq"):
                    fp = tmpdir / f'{record.id}.fq'
                    with fp.open('w') as outfile:
                        SeqIO.write([record], outfile, 'fastq')
                    fpath.append(fp)

        return fpath

    @staticmethod
    def file_len(fname):
        i = 0
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    @staticmethod
    def _diff(res1, res2):
        """ Equal length strings """
        return sum(1 for x, y in zip(res1, res2) if x != y)

    @staticmethod
    def _format_res_string(rstring: str):

        pretty_rstring = ''
        for r in rstring:
            if r.lower() == 'r':
                pretty_rstring += f'{LR}R'
            else:
                pretty_rstring += f'{LB}S'

        return pretty_rstring + f'{Fore.RESET}'

    @staticmethod
    def _format_score(pstring: str):

        try:
            pfloat = float(pstring)
        except ValueError:
            return pstring

        if pfloat < 0.4:
            col = f'{R}'
        elif 0.4 <= pfloat < 0.6:
            col = f'{Y}'
        else:
            col = f'{G}'

        return col + f'{pfloat:.5f}' + f'{Fore.RESET}'

    @staticmethod
    def _natural_key(string_):
        """See http://www.codinghorror.com/blog/archives/001018.html"""
        return [
            int(s) if s.isdigit() else s
            for s in re.split(r'(\d+)', string_)
        ]
