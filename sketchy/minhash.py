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
from uuid import uuid4
from tqdm import tqdm
from pathlib import Path
import dateutil.parser
import multiprocessing as mp

from collections import Counter

from sketchy.utils import PoreLogger, run_cmd

from io import StringIO
import re

from colorama import Fore

LG = Fore.LIGHTGREEN_EX
G = Fore.GREEN
Y = Fore.YELLOW
LR = Fore.LIGHTRED_EX
RE = Fore.RESET
LC = Fore.LIGHTCYAN_EX
C = Fore.CYAN


class MashScore(PoreLogger):

    """ MinHash (MASH) scores for online nanopore predictions """

    def __init__(self):

        PoreLogger.__init__(self)

        # Lachlan suggested keeping a sum od total shared hashes;
        # this which is the default from Sketchy
        self.inter = pandas.DataFrame()  #  MASH outputs, updated
        self.interim = pandas.DataFrame()  # MASH outputs, updated, with data

        self.lineage = Counter()        # Prime lineage counter
        self.genotype = dict()          # Genotype counters by lineage
        self.susceptibility = dict()    # Susceptibility counters by lineage

        self.continuous = list()
        self.start_time_regex = r'start_time=(.*)Z'

    def run(
        self,
        fastq: Path,
        sketch: Path,
        data: Path,
        nreads: int = 100,
        sep: str = '\t',
        index: int = 0,
        cores: int = 2,
        ncpu: int = 0,
        top: int = 10,
        out: Path = None,
        mode: str = 'single',
        sort_by: str = 'shared',
        tmpdir: Path = Path().cwd() / 'tmp',
        show_top: int = 5,
        show_genotype: bool = False,
        nextflow: bool = False,
    ) -> pandas.DataFrame or list:

        df = MashSketch().read_data(fpath=data, sep=sep, index=index)

        if mode == "direct":
            mash = self.mash_dist(
                fastq, mashdb=sketch, ncpu=cores, sort_by=sort_by
            )

            # TODO: Warning: setting index to ID, requires IDs to be UUID

            result_df = mash.drop(
                columns=['file']
            ).set_index('id')[:top]

            result_df.index = result_df.index.map(
                lambda x: Path(x).stem
            )

            print(
                result_df.join(df, how='inner').to_string(header=False)
            )

        else:

            if nextflow:
                # Cut and predict on a single read for dist compute;
                # given by --reads argument on CLI, used in Nextflow
                fpath = self.cut_read(
                    nread=nreads, mode=mode, fastq=fastq, tmpdir=tmpdir
                )
                self.compute_ssh(
                    fpath,
                    df,
                    sketch,
                    nreads,
                    cores,
                    mode,
                    sort_by,
                    tmpdir,
                    show_top,
                    show_genotype,
                    sequential=False
                )

                # Clean up temporary read file
                fpath.unlink()

                return

            if ncpu > 0:

                with mp.Pool(processes=ncpu) as pool:
                    results = [pool.apply_async(
                        func=self._multi_compute, args=(
                            fastq,
                            df,
                            sketch,
                            nread,
                            cores,
                            mode,
                            sort_by,
                            tmpdir,
                            show_top,
                            show_genotype,
                        )) for nread in range(nreads)
                    ]
                    output = [p.get() for p in results]

                    return output
            else:
                for nread in range(nreads):
                    fpath = self.cut_read(
                        nread=nread, mode=mode, fastq=fastq, tmpdir=tmpdir
                    )
                    self.compute_ssh(
                        fpath,
                        df,
                        sketch,
                        nread,
                        cores,
                        mode,
                        sort_by,
                        tmpdir,
                        show_top,
                        show_genotype,
                        sequential=True
                    )

                    # Clean up temporary read file
                    fpath.unlink()
            if out:
                self.inter.to_csv(out, index=True, header=True, sep='\t')

            return self.inter

    def _multi_compute(
            self,
            fastq,
            df,
            sketch,
            nread,
            cores,
            mode,
            sort_by,
            tmpdir,
            show_top,
            show_genotype,
    ):
        """ Wrapper for multiprocessing sum of shared hashes """

        fpath = self.cut_read(
            nread=nread, mode=mode, fastq=fastq, tmpdir=tmpdir
        )
        self.compute_ssh(
            fpath,
            df,
            sketch,
            nread,
            cores,
            mode,
            sort_by,
            tmpdir,
            show_top,
            show_genotype,
            sequential=False
        )

        # Clean up temporary read file
        fpath.unlink()

        return nread

    @staticmethod
    def cut_read(nread, mode, fastq, tmpdir):

        n = 4 * (nread + 1)

        if mode == "cumulative":
            fpath = tmpdir / f'reads_{n}.fq'
            run_cmd(
                f'head -n {n} {fastq} > {fpath}', shell=True
            )
        elif mode == "single":
            fpath = tmpdir / f'read_{n // 4}.fq'
            run_cmd(
                f'head -n {n} {fastq} | tail -4 > {fpath}', shell=True
            )
        else:
            raise

        return fpath

    def compute_ssh(
        self,
        fpath: Path,
        df: pandas.DataFrame,
        sketch: Path,
        nread: int,
        cores: int,
        mode: str = 'single',
        sort_by: str = 'shared',
        tmpdir: Path = Path().cwd() / 'tmp',
        show_top: int = 5,
        show_genotype: bool = False,
        sequential: bool = True,
    ):

        mash = self.mash_dist(
            fpath, mashdb=sketch, ncpu=cores, sort_by=sort_by
        )

        # TODO: NEW METHOD START HERE
        # TODO: Warning: setting index to ID, requires IDs to be UUID

        inter = mash.drop(
            columns=['file', 'dist', 'p-value']
        ).set_index('id')

        if sequential:
            # Update scores in sequential mode
            if self.inter.empty:
                self.inter = inter
            else:
                # Sum of shared hashes
                self.inter = self.inter.add(inter).sort_values(
                    by='shared', ascending=False
                )
            interim = self.inter.copy()
        else:
            # Do not update in distributed mode, output only to file
            interim = inter.copy()

        # Ops to merge with meta data for printing to console
        interim.index = interim.index.map(
            lambda x: Path(x).stem
        )
        interim.index.name = 'uuid'
        interim = interim.join(df, how='inner')

        self.pretty_print(
            interim, fpath, mode, nread, show_top, show_genotype
        )

        # Output sum of shared hashes (
        n = 4 * (nread + 1)
        read = n if mode == "cumulative" else n // 4

        output_prefix = 'total.counts' if sequential else 'read'

        interim.to_csv(
            tmpdir / f'{output_prefix}.{read}', sep='\t'
        )

        if sequential:
            self.interim = interim

        # Legacy method depracated: Select best results from mash dist
        # ordered by shared hash matches, vulnerable to strain mixtures

        # sketchy = None
        # if sketchy == "legacy":
        #     # Originally used with top 2
        #     top_results = mash[:top]
        #
        #     scores = True
        #     if scores:
        #         row = self._compute(
        #             i=nread,
        #             data=df,
        #             tops=top_results,
        #             read_file=fpath,
        #         )
        #         results.append(row)
        #     else:
        #         idx = Path(
        #             top_results.id.values[0]
        #         ).stem
        #         result = df.loc[idx, :]
        #         results.append(result)

        return tmpdir / f'total.counts.{read}'


    def pretty_print(
        self,
        interim: pandas.DataFrame,
        fpath: Path,
        mode: str,
        nread: int,
        select_top: int = 5,
        show_genotype: bool = False,
    ):

        seqlen, timestamp = self._parse_read_stats(fpath)

        lineage_string = f"Lineage match: {LG}" \
                         f"{interim[:1].lineage.tolist()[0]}{G} : " + \
                         f"{' : '.join([str(lineage) for lineage in interim[1:select_top].lineage.tolist()])}"

        suscept_string = f"{Y}Suceptibility: " \
            f"{LC}{interim[:1].susceptibility.tolist()[0]}{RE}"

        geno_string = f"{Y}Genotype: " \
            f"{LC}{interim[:1].genotype.tolist()[0]}{RE}"

        info_string = f"{Y}SSH ({LR}{mode}{Y}) " \
            f"@ read {LR}{str(nread)}{Y}"

        length_string = f"{Y}Read length: {C}{seqlen}"
        time_string = f"{Y}Time: {C}{timestamp}{RE}"

        print(
            f"{info_string:<55}"
            f"{lineage_string:<100}" +
            f"{geno_string if show_genotype else suscept_string:<50}" +
            f"{length_string:<30}" +
            f"{time_string:<30}"
        )
    @staticmethod
    def mash_dist(file, mashdb, ncpu=4, sort_by='shared'):

        result = run_cmd(
            f'mash dist -p {ncpu} {mashdb} {file}', shell=True
        )

        df = pandas.read_csv(
            StringIO(result.decode("utf-8")), sep='\t', header=None,
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
            raise ValueError(
                'MASH distance must be sorted by one of: shared, dist'
            )

    def _compute(
        self,
        i: int,
        data: pandas.DataFrame,
        tops: pandas.DataFrame,
        weight: float = 0.1,
        read_file: Path or str = None,
        quiet: bool = False,
        raw: bool = False,
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

            if raw:
                print(
                    tops.to_string(header=False)
                )
            else:
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

    def _parse_read_stats(self, read_file):

        # last read
        read = run_cmd(
            f'tail -n 4 {read_file}', shell=True
        )

        lines = read.decode("utf-8").split('\n')
        header = lines[0]
        seq = lines[1]

        try:
            timestr = re.search(self.start_time_regex, header)
            if timestr:
                time = timestr.group(1).strip().replace('start_time=', '')
                dtime = dateutil.parser.parse(time)
            else:
                dtime = '-'
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
                pretty_rstring += f'{LC}S'

        return pretty_rstring + f'{Fore.RESET}'

    @staticmethod
    def _format_score(pstring: float):

        return f'{pstring:.5f}'


class MashSketch(PoreLogger):

    """ MinHash (MASH) database sketches for sequence collections """

    def __init__(self):

        PoreLogger.__init__(self)

        self.data: pandas.DataFrame = pandas.DataFrame(None)

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

    def create_refdb_sketch(self, refdb: Path, outdir: Path):
        """ Create a minimal data table from a recursive directory tree
        of the NanoPath RefDB module to quickly create a sketch from
        all *.fasta files in the domain directories directory
        or directory tree representing RefSeq and other reference data
        for taxonomic identification

        This will parse the filenames as:
            <species_taxid>_<lineage_taxid>_<identifier>.fasta

        This can then be used to create a sketch with MASH to run
        MASH dist against for quick taxnomic identification of uncorrected
        nanopore reads.

        """

        pass

    def link(
        self,
        fdir: Path,
        rename: dict = None,
        symlink: bool = True,
        progbar: bool = True,
        uuid: bool = True,
        uuid_file: Path = Path('sketchy.uuid.tsv')
    ) -> [Path]:

        """ Symlink FASTA files into a directory for sketching """

        if fdir.exists():
            print('File path exists, globbing files:', fdir)
            files = list(
                fdir.glob('*')
            )
            uuids = [f.stem for f in files]
        else:
            fdir.mkdir(parents=True)

            if uuid:
                rename = {index: str(uuid4()) for index in self.data.index}

            files = []
            uuids = []
            for i, index in tqdm(
                enumerate(self.data.index),
                disable=not progbar,
                total=len(self.data)
            ):
                fasta = self.data.fasta[i]

                if rename:
                    try:
                        name = rename[index]
                    except TypeError:
                        self.logger.error(
                            f'Could not find {index} in '
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
                uuids.append(name)

        if uuid:
            data = self.data.copy()
            data['id'] = data.index
            data['uuid'] = uuids

            data = data.set_index('uuid')
            data.to_csv(uuid_file, sep='\t')

            print(data)

            key = data.loc[:, ['fasta', 'id']]

            values = data.loc[:, [
                c for c in data.columns if c not in ('fasta', 'id')
            ]]
            key.to_csv(uuid_file.with_suffix('.key'), sep='\t')
            values.to_csv(uuid_file.with_suffix('.data'), sep='\t')

        return files

    @staticmethod
    def sketch(
        fdir: Path,
        name: Path = 'sketchy',
        k: int = 15,
        size: int = 1000,
        glob: str = "*.fasta",
        ncpu: int = 4
    ) -> Path:

        """ Sketch a collection of FASTA files """

        run_cmd(
            f'mash sketch -p {ncpu} -s {size}'
            f' -k {k} -o {name} {fdir}{os.sep}{glob}',
            shell=True
        )

        return name.with_suffix('.msh')

    @staticmethod
    def _check_dir(fdir: str or Path) -> Path:

        fdir = Path(fdir) if isinstance(fdir, str) else fdir

        if not fdir.exists():
            fdir.mkdir(parents=True)

        return fdir





