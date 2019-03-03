"""
=========================================
Module for database sketch using MinHash
=========================================

Prototype for Staphylococcus aureus, Mycobacterium
tuberculosis, Klebsiella pneumoniae.

"""

import os
import shutil
import pandas
import delegator

from tqdm import tqdm
from pathlib import Path

from sketchy.utils import PoreLogger


class MashSketch(PoreLogger):

    """ MinHash (MASH) database sketches for sequence collections """

    def __init__(self,):

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

                Susceptibility pattern from Mykrobe, concatenated into
                a string of length 12, e.g. SSSSSSSSSSSS.

            .. fasta, str:

                Path of FASTA file containing the genome assembly.

        """

        self.data = pandas.read_csv(
            fpath, sep=sep, index_col=index, header=0
        )

    def link(
        self,
        fdir: str or Path,
        rename: dict = None,
        symlink: bool = True,
        progbar: bool = True,
    ):

        """ Symlink FASTA files into a directory for sketching """

        fdir = self._check_dir(fdir)

        for fasta in tqdm(
                self.data.fasta,
                disable=not progbar,
                total=len(self.data)
        ):
            if rename:
                try:
                    name = rename[Path(fasta).stem] + '.fasta'
                except TypeError:
                    self.logger.error(
                        f'Could not find {Path(fasta).stem} in '
                        f'dictionary keys for renaming sketch files.'
                    )
                    continue
            else:
                name = Path(fasta).name

            if symlink:
                (fdir / name).symlink_to(fasta)
            else:
                shutil.copy(
                    fasta, str(fdir / name)
                )

    @staticmethod
    def sketch(
            name: str or Path = 'sketchy',
            fdir: str or Path = Path.cwd(),
            k: int = 15,
            glob="*.fasta"
    ) -> Path:

        """ Sketch a collection of FASTA files """

        name, fdir = Path(name).resolve(), Path(fdir).resolve()

        delegator.run(
            f'mash sketch -k {k} -o {name} {fdir}{os.sep}{glob}'
        )

        return name.with_suffix('.msh')

    @staticmethod
    def _check_dir(fdir: str or Path) -> Path:

        fdir = Path(fdir) if isinstance(fdir, str) else fdir

        if not fdir.exists():
            fdir.mkdir(parents=True)

        return fdir





