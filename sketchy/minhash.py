"""
=========================================
Module for database sketch using MinHash
=========================================

Prototype for Staphylococcus aureus, Klebsiella pneumoniae,
Mycobacterium tuberculosis

"""

import os
import shutil
import pandas
import datetime
import delegator

from tqdm import tqdm
from pathlib import Path

from collections import Counter

from sketchy.utils import PoreLogger

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


class MashScore(PoreLogger):

    """ MinHash (MASH) scores for online nanopore predictions """

    def __init__(self):

        PoreLogger.__init__(self)

        self.lineage = Counter()        # Prime lineage counter
        self.genotype = dict()          # Genotype counters by lineage
        self.susceptibility = dict()    # Susceptibility counters by lineage

        self.continuous = list()

    def run(
        self,
        read_files: [],
        sketch: Path,
        data: Path,
        sep: str = '\t',
        index: int = 0,
        cores: int = 8,
        top: int = 1,
        out: Path = None,
        score: bool = True,
        sort_by: str = 'shared',
        quiet: bool = False,
    ) -> pandas.DataFrame:

        df = MashSketch().read_data(fpath=data, sep=sep, index=index)

        results = []
        for i, read_file in enumerate(read_files):

            mash = self.mash_dist(
                read_file, mashdb=sketch, ncpu=cores, sort_by=sort_by
            )

            # Select best results from mash dist
            # ordered by shared hash matches:

            # TODO: Skip if no shared, maybe something we can do with p-values?
            top_results = mash[:top]

            if score:
                row = self._compute(
                    i=i, data=df, tops=top_results, read_file=read_file, quiet=quiet
                )
                results.append(row)
            else:
                results.append(top_results)

        df = pandas.concat(results, axis=1).T

        if not score:
            df['id'] = df['id'].apply(lambda x: Path(x).stem)
            df['read'] = df['file'].apply(lambda x: Path(x).stem)
            df = df.drop('file', axis=1)

        if out:
            df.to_csv(out, index=False, header=True)

        return df

    @staticmethod
    def mash_dist(file, mashdb, ncpu=4, sort_by='shared'):

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
        df.dist = df.dist.astype(float)

        if sort_by == 'shared':
            return df.sort_values(by=sort_by, ascending=False)
        elif sort_by == 'dist':
            return df.sort_values(by=sort_by, ascending=True)
        else:
            raise ValueError('MASH distance must be sorted by one of: shared, dist')

    def _compute(
            self,
            i: int,
            data: pandas.DataFrame,
            tops: pandas.DataFrame,
            weight: float = 0.1,
            read_file: Path or str = None,
            quiet: bool = False
    ) -> pandas.Series:
        """ Update counts and compute preference score

        :param i: index of read, number of read in cumulative file
        :param tops: dataframe of top shared hashes against reference sketch
        :param weight: modifier for consecutive lineage weight
        :param read_file: file containing reads, for parsing last read stats
        :return:
        """

        # Update observed ST and resistance profile within ST
        for tid in tops.id:
            iid = Path(tid).stem
            st = str(data.at[iid, 'lineage'])
            gen = str(data.at[iid, 'genotype'])
            res = str(data.at[iid, 'susceptibility'])

            self.lineage.update([st])

            if st not in self.genotype.keys():
                self.genotype[st] = Counter()
            if st not in self.susceptibility.keys():
                self.susceptibility[st] = Counter()

            self.genotype[st].update([gen])
            self.susceptibility[st].update([res])

        # Top 2 lineages for score:
        top_lineages = self.lineage.most_common(2)   # in all STs

        top_st = top_lineages[0][0]
        top_count = top_lineages[0][1]

        try:
            second_st = top_lineages[1][0]
            second_count = top_lineages[1][1]
        except IndexError:
            second_st, second_count = "", ""

        # Top within lineage susceptibility profile:
        top_within_lineage_susceptibility = \
            self.susceptibility[top_st].most_common(1)[0][0]  # in top lineage

        # Top within lineage genotype profile:
        top_within_lineage_genotype = \
            self.genotype[top_st].most_common(1)[0][0]  # in top lineage

        # Weighted score based on calling the
        # same lineage consecutively:
        if not self.continuous:
            self.continuous.append(top_st)
        else:
            if self.continuous[-1] == top_st:
                self.continuous.append(top_st)
            else:
                self.continuous = list()

        # PSG like preference score, see Brinda et al. 2019
        # Weighted score can force a preference score > 1
        # so the preference score was bounded the score at 1:
        try:
            score = 2*top_count/(second_count + top_count) - 1

            if weight:
                score += (
                    weight * len(self.continuous)
                )

            if score > 1.:
                score = 1.

        except TypeError:
            # If second count not available on first prediction,
            # i.e. when second count raises IndexError above:
            score = 0.

        seqlen, time = self._parse_read_stats(read_file)

        if not quiet:
            print(
                f"{i}\t",
                f"{top_st}\t",
                f"{top_count}\t",
                f"{second_st}\t",
                f"{second_count}\t",
                f"{self._format_score(score)}\t",
                f"{top_within_lineage_susceptibility}\t",
                f"{top_within_lineage_genotype}\t"
                f"{seqlen}\t",
                f"{time}"
            )

        return pandas.Series(
            data={
                'primary_lineage': top_st,
                'primary_count': top_count,
                'secondary_lineage': second_st,
                'secondary_count': second_count,
                'score': self._format_score(score),
                'suscptibility': top_within_lineage_susceptibility,
                'genotype': top_within_lineage_genotype,
                'last_read_length': seqlen,
                'last_read_start_time': time
            }
        )

    @staticmethod
    def _parse_read_stats(read_file):

        # last read
        read = delegator.run(
            f'tail -n 4 {read_file}'
        )

        lines = read.out.split('\n')
        header = lines[0]
        seq = lines[1]

        try:
            time = header.split('start_time=')[1]
            time = time.replace('T', '-').strip('Z')

            dtime = datetime.datetime.strptime(time, '%Y-%m-%d-%H:%M:%S')
        except IndexError:
            dtime = '-'

        return len(seq), str(dtime)

    @staticmethod
    def _format_res_string(rstring: str):

        if rstring in ('nan', 'None'):
            return '-'

        pretty_rstring = ''
        for r in rstring:
            if r.lower() == 'r':
                pretty_rstring += f'{LR}R'
            else:
                pretty_rstring += f'{LB}S'

        return pretty_rstring + f'{Fore.RESET}'

    @staticmethod
    def _format_score(pstring: float):

        return f'{pstring:.5f}'


class MashSketch(PoreLogger):

    """ MinHash (MASH) database sketches for sequence collections """

    def __init__(self):

        PoreLogger.__init__(self)

        self.data: pandas.DataFrame = None

    def read_data(self, fpath, sep='\t', index=0):

        """ Read a data frame with all information necessary to
        sketch the genomes with MASH. This includes one row per
        genome over the following column headers, including an
        index of isolate identifiers:

            .. index, str:

                Data frame index column for isolate identifiers,
                e.g. ERR2547368 for an isolate downloaded from
                the European Nucleotide Archive

            .. lineage, str:

                Lineage classification for the organism of interest,
                for instance MLST in Staphylococcus aureus or SNP
                types for Mycobacterium tuberculosis.

            .. susceptibility, str:

                Antibiotic resistance and susceptibility patterns, e.g. from
                Mykrobe concatenated into a string.

            .. genotype, str:

                Genotype, e.g. typing pattern from Kleborate,
                concatenated into a string.

            .. fasta, str:

                Path of FASTA file containing the genome assembly.

        """

        self.data = pandas.read_csv(
            fpath, sep=sep, index_col=index, header=0
        )

        return self.data

    def link(
        self,
        fdir: str or Path,
        rename: dict = None,
        symlink: bool = True,
        progbar: bool = True,
    ) -> [Path]:

        """ Symlink FASTA files into a directory for sketching """

        fdir = self._check_dir(fdir)

        files = []
        for fasta in tqdm(
                self.data.fasta,
                disable=not progbar,
                total=len(self.data)
        ):
            if rename:
                try:
                    name = rename[Path(fasta).stem]
                except TypeError:
                    self.logger.error(
                        f'Could not find {Path(fasta).stem} in '
                        f'dictionary keys for renaming sketch files.'
                    )
                    continue
            else:
                name = Path(fasta).stem

            if symlink:
                (fdir / name).symlink_to(fasta)
            else:
                shutil.copy(
                    fasta, str(fdir / name)
                )

            files.append(fdir / name)

        return files

    @staticmethod
    def sketch(
        name: str or Path = 'sketchy',
        fdir: str or Path = Path.cwd(),
        k: int = 15,
        size: int = 1000,
        glob="*.fasta"
    ) -> Path:

        """ Sketch a collection of FASTA files """

        name, fdir = Path(name).resolve(), Path(fdir).resolve()

        delegator.run(
            f'mash sketch -s {size} -k {k} -o {name} {fdir}{os.sep}{glob}'
        )

        return name.with_suffix('.msh')

    @staticmethod
    def _check_dir(fdir: str or Path) -> Path:

        fdir = Path(fdir) if isinstance(fdir, str) else fdir

        if not fdir.exists():
            fdir.mkdir(parents=True)

        return fdir





