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

        self.inter = pandas.DataFrame()    # MASH outputs, updated
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
        mode: str = 'single',
        tmpdir: Path = Path().cwd() / 'tmp',
        show_top: int = 5,
        show_genotype: bool = False,
        pretty: bool = False,
        info: bool = False,
    ) -> pandas.DataFrame or list:

        df = MashSketch().read_data(fpath=data, sep=sep, index=index)

        total_reads = self._get_total_reads(fastq)

        if mode == "direct":
            # Direct mode on file:
            mash = self.mash_dist(
                fastq, mashdb=sketch, ncpu=cores, sort_by='shared'
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
            # Multiprocessing file mode:
            if ncpu > 0:
                if nreads > total_reads:
                    if nreads > total_reads:
                        self.logger.info(
                            f'Selected reads {nreads} are larger than'
                            f' total number of reads: {total_reads}'
                        )
                        self.logger.info(
                            f'Set reads to total number of reads: {total_reads}'
                        )
                        nreads = total_reads

                with mp.Pool(processes=ncpu) as pool:
                    results = [pool.apply_async(
                        func=self._multi_compute, args=(
                            fastq,
                            df,
                            sketch,
                            nread,  # print read at >> nread << in loop
                            cores,
                            mode,
                            tmpdir,
                        ), callback=lambda x: self.logger.debug(
                            f'Computed shared hashes for read '
                            f'{x.stem.split("_")[-1]} @ {x}'
                        ))
                        for nread in range(nreads)
                    ]

                    output = [p.get() for p in tqdm(results)]

                    return output
            else:
                # Online mode, single CPU
                for nread in range(nreads):
                    if nread == total_reads:
                        self.logger.info(
                            f'Reached total number of '
                            f'reads in file: {total_reads}'
                        )
                        # Exit gracefully prevents braking Nextflows
                        exit(0)

                    fpath = self.cut_read(
                        nread=nread, mode=mode, fastq=fastq, tmpdir=tmpdir
                    )
                    self.compute_ssh(
                        fpath=fpath,
                        df=df,
                        sketch=sketch,
                        nread=nread,  # print read at >> nread << in loop
                        cores=cores,
                        mode=mode,
                        sort_by='shared',
                        tmpdir=tmpdir,
                        show_top=show_top,
                        show_genotype=show_genotype,
                        sequential=True,
                        pretty=pretty,
                        info=info,
                    )

                    # Clean up temporary read file
                    fpath.unlink()

                    return self.inter


    def _multi_compute(
            self,
            fastq,
            df,
            sketch,
            nread,
            cores,
            mode,
            tmpdir,
    ):
        """ Wrapper for multiprocessing sum of shared hashes """

        fpath = self.cut_read(
            nread=nread, mode=mode, fastq=fastq, tmpdir=tmpdir
        )
        self.compute_ssh(
            fpath=fpath,
            df=df,
            sketch=sketch,
            nread=nread,
            cores=cores,
            mode=mode,
            sort_by='shared',
            tmpdir=tmpdir,
            sequential=False,
            pretty=False,
            info=False,
            printout=False,
        )

        fpath.unlink()

        return fpath

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

    def _get_total_reads(self, fastq: Path) -> int or None:

        try:
            out = run_cmd(f'wc -l {fastq}', shell=True)
            return int(out.decode("utf-8").strip('\n').split()[0]) // 4
        except:
            self.logger.debug('Error in getting total read count.')
            self.logger.info(f'Could not detect total read number in {fastq}')
            raise

    def compute_ssh(
        self,
        fpath: Path = None,
        sketch: Path = None,
        nread: int = None,
        cores: int = 2,
        mash: pandas.DataFrame or None = None,
        df: pandas.DataFrame or None = None,
        mode: str = 'single',
        sort_by: str = 'shared',
        tmpdir: Path = Path().cwd() / 'tmp',
        show_top: int = 5,
        show_genotype: bool = False,
        sequential: bool = True,
        pretty: bool = False,
        info: bool = False,
        printout: bool = True,
    ):

        if mash is None:
            mash = self.mash_dist(
                fpath, mashdb=sketch, ncpu=cores, sort_by=sort_by
            )

            inter = mash.drop(
                columns=['file', 'dist', 'p-value']
            ).set_index('id')
        else:

            inter = mash.drop(
                columns=['lineage', 'genotype', 'susceptibility']
            )  # not setting index done before on mash: .set_index('uuid')

        if sequential:
            # Update scores in sequential mode
            if self.inter.empty:
                self.inter = inter
            else:
                # Sum of shared hashes - add (original)
                self.inter = self.inter.add(inter, fill_value=0).sort_values(
                    by='shared', ascending=False
                )
            interim = self.inter.copy()
        else:
            # Do not update in distributed mode,
            # output reduced (filtered) to files
            interim = inter.copy()
            interim = interim[interim.shared > 0].sort_values(
                by='shared', ascending=False
            )


        # Ops to merge with meta data for printing to console
        interim.index = interim.index.map(
            lambda x: Path(x).stem
        )
        interim.index.name = 'uuid'

        if df is not None:
            interim = interim.join(df, how='inner')

        if nread is not None:
            # Output sum of shared hashes (
            n = 4 * (nread + 1)
            read = n if mode == "cumulative" else n // 4

            output_prefix = 'total.counts' if sequential else 'read'

            if not interim.empty:
                interim.to_csv(
                    tmpdir / f'{output_prefix}.{read}', sep='\t'
                )

        if printout and not interim.empty:

            if info:
                seqlen, timestamp = self._parse_read_stats(fpath)
            else:
                seqlen, timestamp = None, None

            if pretty:
                self.pretty_print(
                    interim=interim,
                    mode=mode,
                    nread=nread,
                    seqlen=seqlen,
                    timestamp=timestamp,
                    select_top=show_top,
                    show_genotype=show_genotype
                )
            else:
                self.regular_print(
                    interim=interim,
                    nread=nread,
                    seqlen=seqlen,
                    timestamp=timestamp,
                    select_top=show_top,
                    show_genotype=show_genotype
            )

        if sequential:
            self.interim = interim

    @staticmethod
    def regular_print(
        interim: pandas.DataFrame,
        nread: int,
        seqlen: int,
        timestamp: str,
        select_top: int = 5,
        show_genotype: bool = False,
    ):

        lineages = ':'.join([
            str(lineage) for lineage in interim[:select_top].lineage.tolist()
        ])
        susceptibility = interim[:1].susceptibility.tolist()[0]
        genotype = interim[:1].genotype.tolist()[0]

        print(
            f"{nread}\t"
            f"{lineages}\t"
            f"{seqlen}\t",
            f"{timestamp}\t",
            f"{susceptibility}\t",
            f"{genotype if show_genotype else None}"
        )

    @staticmethod
    def pretty_print(
        interim: pandas.DataFrame,
        mode: str,
        nread: int,
        seqlen: int,
        timestamp: str,
        select_top: int = 5,
        show_genotype: bool = False,
    ):

        lineages = ' : '.join([
            str(lineage) for lineage in interim[1:select_top].lineage.tolist()
        ])

        lineage_string = f"Lineage match: {LG}" \
            f"{interim[:1].lineage.tolist()[0]}{G} : " + \
            f"{lineages}"

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

        """ Compute MASH distance by simply calling MASH

        :param file:
        :param mashdb:
        :param ncpu:
        :param sort_by:
        :return:
        """

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

    def _parse_read_stats(self, read_file):

        # last read
        read = run_cmd(
            f'tail -n 4 {read_file}', shell=True
        )

        lines = read.decode("utf-8").split('\n')

        try:
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

        except KeyError:
            self.logger.info(f'Could not detect last read in: {read_file}')
            dtime, seq = '-', ''

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





