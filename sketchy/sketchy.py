"""
================================
Access point module for Sketchy
================================

TODO: do not jointly count genotypes, count each individually!

"""

import re
import pandas
import random
import numpy

import delegator

from pathlib import Path

from sketchy.minhash import MashScore
from sketchy.minhash import MashSketch

from Bio import SeqIO
from datetime import datetime
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

    @staticmethod
    def create_mash_sketch(
        data: Path or str = Path().home() / 'data.tab',
        kmer_length: int or [int] = 15,
        sketch_size: int or [int] = 1000,
        outdir: Path or str = Path().home(),
        file_glob: str = "*.fasta",
        file_copy: bool = False,
        prefix: str = 'sketchy',
        sep: str = '\t',
        index_col: int = 0
    ):

        sketch = MashSketch()

        sketch.read_data(
            fpath=data, sep=sep, index=index_col
        )

        files = sketch.link(
            fdir=outdir,
            rename=None,
            symlink=not file_copy,
            progbar=True
        )

        sketch.sketch(
            name=f'{prefix}_{kmer_length}',
            fdir=outdir,
            k=kmer_length,
            size=sketch_size,
            glob=file_glob
        )

        for file in files:
            file.unlink()

    def predict_nanopore(
            self,
            fastq: Path or str,
            sketch: Path or str,
            data: Path or str,
            cores: int = 16,
            score: bool = True,
            header: bool = False,
            nreads: int = 100,
            top: int = 1,
            out: Path = Path().cwd() / 'results.csv',
            tmp: Path = Path.cwd() / 'tmp',
            sort_by: str = 'shared'
    ) -> pandas.DataFrame:

        """ Online predictions on nanopore reads

        :param fastq:
            fastq file to slice into temporary cumulative slices to
            compute scores on and simulate sequencing, extracts
            timestamps from read header (Albacore)

        :param sketch:
            mash sketch of database to compute prediction on

        :param cores:
            number of compute cores for mash dist

        :param top_results:
            number of top mash dist results to compute scores with

        :returns
            None
        """

        tmp.mkdir(parents=True, exist_ok=True)

        if score:
            # Cumulative read slicing for score predictions
            reads = [
                self._slice_fastq(
                    fastq=fastq, nreads=i, tmpdir=tmp, cumulative=True
                ) for i in range(1, nreads+1, 1)
            ]
        else:
            # Single read slicing for shared hashes output:
            reads = self._slice_fastq(
                fastq=fastq, nreads=1, tmpdir=tmp, cumulative=False
            )

        if header:
            self._print_header1()

        ms = MashScore()
        results = ms.run(
            read_files=reads,
            sketch=sketch,
            cores=cores,
            top=top,
            score=score,
            data=data,
            out=out,
            sort_by=sort_by
        )

        return results

    @staticmethod
    def _print_header1():

        print(
            f"{C}{'-' * 143}{RE}\n"
            f"{LC}{'Read':<5}{RE}",
            f"{LM}{'ST:1':<7}{RE}",
            f"{LC}{'Count':<7}{RE}",
            f"{LM}{'ST:2':<7}{RE}",
            f"{LC}{'Count':<7}{RE}",
            f"{LC}{'Score':<7}{RE}",
            f"{LY}{'Profile':<15}{RE}",
            f"{LY}{'Genotype':<50}{RE}",
            f"{LC}{'Length':<10}{RE}",
            f"{LC}{'Start Time':<15}{RE}",
            f"\n{C}{'-' * 143}{RE}"
        )

    @staticmethod
    def _print_header2():

        print(
            f"{C}{'-' * 75}{RE}\n"
            f"{LC}{'Genome':<7}{RE}",
            f"{LM}{'Predict':<10}{RE}",
            f"{LY}{'Predict':<12}{RE}",
            f"{LC}{'Hashes':<10}{RE}",
            f"{LM}{'True':<10}{RE}",
            f"{LY}{'True':<12}{RE}",
            f"{LY}{'Diff':<5}{RE}",
            f"\n{C}{'-' * 75}{RE}"
        )

    @staticmethod
    def sort_fastq(file: str = None, fastq: str = None, shuffle: bool = False, nbootstrap: int = None,
                   sample_size: int = None, replacement: bool = False, prefix: str = 'boot_') -> list or str:
        """ Not very efficient FASTQ sort and shuffle """

        dates = []
        ids = []
        lengths = []
        records = {}
        with open(file, "r") as input_handle:
            for record in SeqIO.parse(input_handle, "fastq"):
                # Extract start time
                try:
                    time = record.description.split('start_time=')[1]
                    time = time.replace('T', '-').strip('Z')
                    dtime = datetime.strptime(time, '%Y-%m-%d-%H:%M:%S')
                except IndexError:
                    time = '-'
                    dtime = '-'

                dates.append(dtime)

                ids.append(record.id)
                lengths.append(
                    len(record.seq)
                )

                records[record.id] = record

                if not fastq:
                    print(record.id, time, len(record.seq))

        df = pandas.DataFrame(
            data={
                'date': dates,
                'read': ids,
                'length': lengths
            }
        ).sort_values(by='date')

        pandas.set_option('display.max_rows', 120)

        df = df.reset_index()

        if fastq:

            recs = [
                records[read] for read in df['read']
            ]

            if shuffle:
                recs = random.shuffle(recs)

            if nbootstrap and sample_size:
                rec_samples = [numpy.random.choice(
                    numpy.array(recs),
                    replace=replacement,
                    size=sample_size
                ).tolist() for _ in range(nbootstrap)]

                fnames = []
                for i, rec in enumerate(rec_samples):
                    fname = prefix + str(i) + '.fq'
                    with open(fname, "w") as output_handle:
                        SeqIO.write(recs, output_handle, 'fastq')
                    fnames.append(fname)

                return fnames
            else:
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
            delegator.run(
                f'head -n {nreads} {fastq.resolve()} > {fpath.resolve()}'
            )
        else:
            fpath = []
            with fastq.resolve().open('r') as input_handle:
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
