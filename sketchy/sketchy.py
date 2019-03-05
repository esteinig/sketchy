"""
================================
Access point module for Sketchy
================================


"""

import re
import pandas
import random

import delegator

import os

from pathlib import Path
from collections import Counter
from io import StringIO

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

from sketchy.minhash import MashSketch


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
            outdir: Path or str = Path().home(),
            file_glob: str = "*.fasta",
            file_copy: bool = False,
            prefix: str = 'sketchy',

    ):

        sketch = MashSketch()

        # Read Data Frame with columns:
        #   IndexID
        sketch.read_data(
            fpath=data, sep='\t', index=0
        )

        sketch.link(
            fdir=outdir,
            rename=None,
            symlink=not file_copy,
            progbar=True
        )

        sketch.sketch(
            name=f'{prefix}_{kmer_length}',
            fdir=Path.cwd(),
            k=kmer_length,
            glob=file_glob
        )

    def predict_assemblies(
            self,
            assemblies: Path or str,
            sketch: str or Path,
            glob: str = '*.fasta',
            header: bool = True,
            cores: int = 16,
            top_results: int = 1,
    ):

        genomes = Path(assemblies).glob(glob)

        if header:
            self._print_header2()
        for i, genome in enumerate(genomes):
            tops = self.mash_dist(genome, sketch, ncpu=cores, top=top_results)
            self._compute_scores(i, tops, genome=genome.name)

    def predict_nanopore(
            self,
            sketch: Path or str,
            fastq: Path or str = None,
            tmp: Path or str = Path.cwd() / 'tmp',
            watch_dir: Path or str = None,
            cores: int = 16,
            top_results: int = 1,
            header: bool = False,
            nreads: int = 100,

    ) -> None:
        """ Online predictions on nanopore reads

        :param sketch:
            mash sketch of database to compute prediction on

        :param fastq:
            fastq file to slice into temporary cumulative slices to
            compute scores on and simulate sequencing, extracts
            timestamps from read header

        :param watch_dir:
            directory path to watch for new fastq files of reads
            or when read is added to fastq

        :param test_dir:
            directory with cumulative fastq files for simulation

        :param extension:
            extension for iterating over cumulative fastq files

        :param cores:
            number of compute cores for mash dist

        :param top_results:
            number of top dist results to compute scores on

        :returns
            None
        """

        if fastq:

            if not nreads % 4 == 0:
                raise ValueError(f'Fastq format not recognized in {fastq}.')

            tmp.mkdir(parents=True, exist_ok=True)

            reads = [
                self._slice_fastq(
                    fastq=fastq, nreads=i, tmpdir=tmp
                ) for i in range(1, nreads+1, 1)
            ]

        else:
            raise ValueError(
                'Either param `test_dir` or `fastq` must be specified.'
            )

        if header:
            self._print_header1()

        lineage = Counter()
        resistance = Counter()
        continuous = list()
        for i, read in enumerate(reads):
            tops = self.mash_dist(
                read, mashdb=sketch, ncpu=cores, top=top_results
            )

            top = self._compute_scores(
                i, tops, lineage=lineage, resistance=resistance,
                continuous=continuous
            )

            if not continuous:
                continuous.append(top)
            else:
                if top == continuous.pop():
                    continuous.append(top)
                else:
                    continuous = list()

            if fastq:
                # Delete tmp slices
                os.remove(
                    str(read)
                )

    @staticmethod
    def mash_dist(file, mashdb, ncpu=4, top=2):

        result = delegator.run(
            f'mash dist -p {ncpu} {mashdb} {file}'
        )
        df = pandas.read_csv(
            StringIO(result.out), sep='\t', header=None,
            names=[
                "id", 'file', 'dist', "p-value", "shared"
            ], index_col=False
        )

        shared = pandas.DataFrame(
            df.shared.str.split('/').tolist(), columns=['shared', 'total']
        )

        df.shared = shared.shared.astype(int)

        df = df.sort_values(by='shared', ascending=False)

        if top:
            df = df[:top]

        return df

    @staticmethod
    def _print_header1():

        print(
            f"{C}{'-' * 60}{RE}\n"
            f"{LC}{'Read':<5}{RE}",
            f"{LM}{'ST:1':<7}{RE}",
            f"{LC}{'Count':<7}{RE}",
            f"{LY}{'Profile':<12}{RE}",
            f"{LM}{'ST:2':<7}{RE}",
            f"{LY}{'Count':<7}{RE}",
            f"{LY}{'Score':<5}{RE}\n",
            f"{C}{'-' * 60}{RE}"
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
            f"{LY}{'Diff':<5}{RE}\n",
            f"{C}{'-' * 75}{RE}"
        )

    def _compute_scores(
            self,
            i: int,
            tops: pandas.DataFrame,
            lineage: Counter = None,
            resistance: Counter = None,
            genome: str = None,
            continuous: list = None,
            weight: float = 0.1
    ):

            iids, sts, resist, mashshare = [], [], [], []
            for tid in tops.id:
                iid, st, res = tid.strip('.fasta').split('_')

                if lineage is not None and resistance is not None:
                    lineage.update([st])
                    resistance.update([res])
                else:
                    iids.append(iid)
                    sts.append(st)
                    resist.append(res)
                    mashshare.append(tops[tops['id'] == tid].shared.values[0])

            if lineage is not None and resistance is not None:
                lin = lineage.most_common(3)
                rest = resistance.most_common(3)

                top_st = lin[0][0]
                top_count = lin[0][1]

                try:
                    second_st = lin[1][0]
                    second_count = lin[1][1]
                except IndexError:
                    second_st, second_count = "", ""

                try:
                    # PSG like count score, see Brinda et al. 2019
                    # Includes optional weight to add from confident sequential
                    # ST predictions - this can force a preference score > 1
                    # so we bounded the score at 1:
                    ratio = 2*top_count/(second_count + top_count) - 1 + \
                            (weight * len(continuous))

                    if ratio > 1.:
                        ratio = 1.

                except TypeError:
                    ratio = ''

                if isinstance(ratio, float):
                    col = f"{RE}" if ratio < 0.6 else f"{G}"
                else:
                    col = f"{RE}"

                print(
                    f"{i:<5}",
                    f"{col}{'ST' + top_st:<7}{RE}",
                    f"{top_count:<7}",
                    f"{self._format_res_string(rest[0][0]):<15}",
                    f"{'ST' + second_st:<7}",
                    f"{second_count:<7}",
                    f"{self._format_score(ratio):<5}"
                )

                return top_st
            else:
                if genome:
                    giid, gst, gres = genome.strip('.fasta').split('_')
                else:
                    giid, gst, gres = '-', '-', '-'

                topst = sts[0]

                diff = self._diff(resist[0], gres)
                print(
                    f"{i:<7}",
                    f"{R if topst != gst else G}{'ST' + topst:<10}{RE}",
                    f"{self._format_res_string(resist[0]):<15}",
                    f"{mashshare[0]:<10}",
                    f"{'ST' + gst:<10}",
                    f"{self._format_res_string(gres):<15}",
                    f"{R if diff > 0 else G}{diff:<7}{RE}",
                )

    @staticmethod
    def sort_fastq(file, fastq: str = None, shuffle: bool = False):

        from Bio import SeqIO
        from datetime import datetime

        dates = []
        ids = []
        lengths = []
        records = {}
        with open(file, "r") as input_handle:
            for record in SeqIO.parse(input_handle, "fastq"):
                time = record.description.split('start_time=')[1]
                time = time.replace('T', '-').strip('Z')
                dtime = datetime.strptime(time, '%Y-%m-%d-%H:%M:%S')
                dates.append(dtime)
                ids.append(record.id)
                lengths.append(len(record.seq))
                records[record.id] = record
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

        # 3 percent of genome with large database

        if fastq:

            recs = [
                records[read] for read in df['read']
            ]

            if shuffle:
                recs = random.shuffle(recs)

            with open(fastq, "w") as output_handle:
                SeqIO.write(recs, output_handle, 'fastq')

    @staticmethod
    def _slice_fastq(
            fastq: Path,
            nreads: int = 1,
            tmpdir: Path = Path().cwd() / 'tmp'
    ) -> Path:

        fpath = tmpdir / f'reads_{nreads}.fq'

        nreads = 4*nreads

        delegator.run(
            f'head -n {nreads} {fastq.resolve()} > {fpath.resolve()}'
        )

        return fpath

    @staticmethod
    def file_len(fname):
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

        pfloat = round(pfloat, 5)

        return col + str(pfloat) + f'{Fore.RESET}'

    @staticmethod
    def _natural_key(string_):
        """See http://www.codinghorror.com/blog/archives/001018.html"""
        return [
            int(s) if s.isdigit() else s
            for s in re.split(r'(\d+)', string_)
        ]
