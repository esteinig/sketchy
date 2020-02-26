import pandas
import copy
import shutil

from tqdm import tqdm
from pathlib import Path
from dataclasses import dataclass
from pandas.core.groupby import DataFrameGroupBy
from sketchy.utils import PoreLogger

#################################
# Data containers for Nextflows #
#################################


@dataclass
class ResultData:
    """ Methods for data containers. """

    iid: pandas.DataFrame = pandas.DataFrame(
        columns=['iid']
    )
    fasta: pandas.DataFrame = pandas.DataFrame(
        columns=['fasta']
    )
    fastq: pandas.DataFrame = pandas.DataFrame(
        columns=['forward', 'reverse']
    )

    def __copy__(self):
        """ Overwritten in subclasses """

    def __iter__(self):

        for attr, df in self.__dict__.items():
            yield attr, df

    def __add__(self, other):

        for attr, data in self:
            data_other = getattr(other, attr)

            merged = pandas.concat(
                [data, data_other], sort=False
            )

            setattr(self, attr, merged)

        self.update_iid()

        return self

    def __getitem__(self, iid):

        return self.iid[iid]

    def isolate(self, *attrs):

        for attr, df in self:
            if attr not in attrs:
                setattr(self, attr, pandas.DataFrame())

    def update_iid(self, iids=None):

        """ Update list of unique IIDs in current data container """

        if not iids:
            iids = pandas.Series(
                [
                    iid for attr, df in self
                    for iid in df.index.unique().tolist()
                ]
            ).unique()

        self.iid = pandas.DataFrame(data={'iid': iids}, index=iids)

    def intersection(self):

        indices = [set(df.index) for attr, df in self if not df.empty]
        intersect = set.intersection(*indices)

        for attr, df in self:
            if not df.empty:
                df = df[~df.index.duplicated()]
                df = df[df.index.isin(intersect)]
                setattr(self, attr, df)

    def remove(self, remove_dict: dict, retain: bool = False):

        for process, df in self:

            # TODO: Patched, needs proper thought.

            try:
                remove_dict[process]
            except KeyError:
                # Skip self.iid, self.fasta, self.fastq
                continue

            if retain:
                removed_df = df[
                    df.index.isin(
                        remove_dict[process]
                    )
                ]
            else:
                removed_df = df.drop(
                    remove_dict[process], errors='ignore'
                )

            setattr(self, process, removed_df)

        self.update_iid()

    def complete(
        self,
        at: list = ('kraken', 'mlst')
    ):

        """ Return data object with the set intersection between all IIDs
        in the analyses results; usually after filtering, can be narrowed
        to only some results by default to Kraken and MLST """

        if at:
            allowed = at
        else:
            allowed = [process for process, _ in self]

        iids = [
            set(df.index.tolist()) for process, df in self
            if len(df) > 0 and process in allowed
        ]

        for li in iids:
            print('Len IID data:', len(li))

        intersect = list(set.intersection(*iids))

        print('Intersect, keep this many:', len(intersect))

        self.remove(
            {process: intersect for process, _ in self},
            retain=True
        )

    def write(self, outdir: str or Path) -> None:

        outdir = Path(outdir) \
            if isinstance(outdir, str) else outdir

        outdir.mkdir(parents=True, exist_ok=True)

        for process, data in self:
            data.to_csv(
                (outdir / f'{process}.csv'),
                index=True, header=True, sep=','
            )

    def read(self, results: str or Path):

        results = Path(results) \
            if isinstance(results, str) else results

        for file in results.glob('*.csv'):
            setattr(self, file.stem, pandas.read_csv(
                file, index_col=0, header=0, sep=','
            ))

        self.update_iid()

    def subset(self, iids, inplace=False):

        """ Subset the Data by a list of isolate IDs.

        If IDs are not in the data frame index, a row of null values
        for this ID is introduced into the data frame.
        """

        if inplace:
            sub = self
        else:
            sub = copy.copy(self)

        for attr, df in self:
            subset = df[df.index.isin(iids)]
            missing = list(set(iids).difference(set(subset.index)))
            if missing:
                subset = pandas.concat(
                    (
                        subset,
                        pandas.DataFrame(index=missing, columns=df.columns)
                     )
                )

            setattr(sub, attr, subset)

        sub.update_iid()

        return sub

    def groupby(self, attr, column, set_index=False, **kwargs) -> [DataFrameGroupBy]:

        """ Map the column values to each corresponding index
        in all other data attributes and group the data frames
        by the column values. This allows e.g. to group all data
        sequence type etc.
        """

        if not hasattr(self, attr):
            raise AttributeError(
                f'ResultData has no attribute: {attr}'
            )

        # Map sequence type to index in each DataFrame as new column
        # and then group every attribute data frame by this column:

        data = getattr(self, attr)
        data_series = data[column]

        for attr, df in self:
            df[column] = df.index.map(data_series)

            if set_index:
                df = df.set_index([df.index, column])
                setattr(self, attr, df)

            groups = df.groupby(by=column, **kwargs)
            yield groups, attr

    def select(
        self, attr, column, min_count=None, sample=None, values=None
    ):

        data = getattr(self, attr)

        # TODO: Both

        if values:
            iids = data[
                data[column].isin(values)
            ].index.unique().tolist()

        elif min_count:
            counts = data[column].value_counts()
            largest = counts[counts > min_count]

            subset = largest.index.to_series()  # Index: unique values counted

            sub = data[
                data[column].isin(subset)
            ]

            iids = sub.index.unique().tolist()

            if sample:

                iids = pandas.concat([
                    group.sample(n=sample) for i, group in sub.groupby(column)
                ]).index

        else:
            iids = data[column].index.unique().tolist()

        return self.subset(iids)

    def link_fasta(self, fdir: str or Path, symlink: bool = True,
                   progbar: bool = True, index: dict = None):

        fdir = self._check_dir(fdir)

        for fasta in tqdm(
            self.fasta.fasta,
            disable=not progbar,
            total=len(self.fasta)
        ):
            if index:
                try:
                    name = index[Path(fasta).stem] + '.fasta'
                except TypeError:
                    print(index[Path(fasta).stem])
                    continue
                except KeyError:
                    print('Could not find entry in FASTA Index.')
                    continue
            else:
                name = Path(fasta).name

            if symlink:
                (fdir / name).symlink_to(fasta)
            else:
                shutil.copy(
                    fasta, str(fdir / name)
                )

    def link_fastq(self, fdir: str or Path, symlink: bool = True,
                   progbar: bool = True):

        fdir = self._check_dir(fdir)

        for fwd, rv in tqdm(
            self.fastq.itertuples(index=False),
            disable=not progbar,
            total=len(self.fastq)
        ):
            for fastq in (fwd, rv):
                if symlink:
                    (fdir / Path(fastq).name).symlink_to(fastq)
                else:
                    shutil.copy(
                        fastq, str(fdir)
                    )

    @staticmethod
    def _check_dir(fdir: str or Path) -> Path:

        fdir = Path(fdir) if isinstance(fdir, str) else fdir

        if not fdir.exists():
            fdir.mkdir(parents=True)

        return fdir

    # File operations

    def _add_fasta(self, fasta_dir: Path, extension: str = '.fasta'):

        """ Add file paths to data frame in attribute: `fasta` """

        # Clear FASTA, weirdly in loops files accumulate
        # even if new instance of this is created

        self.fasta = pandas.DataFrame()

        for fpath in fasta_dir.glob(f'*{extension}'):
            iid = fpath.name.replace(extension, '')
            self.fasta.at[iid, 'fasta'] = fpath

    def _add_fastq(self, fastq_dir: Path, forward_tail: str = '_1',
                   reverse_tail: str = '_2', extension: str = '.fq.gz'):

        self.fastq = pandas.DataFrame()

        for fpath in fastq_dir.glob(f'*{forward_tail}{extension}'):
            iid = fpath.name.replace(f'{forward_tail}{extension}', '')
            if iid in self.fastq.index:
                self.fastq.at[iid, 'forward'] = fpath
            else:
                pass

        for fpath in fastq_dir.glob(f'*{reverse_tail}{extension}'):
            iid = fpath.name.replace(f'{reverse_tail}{extension}', '')
            if iid in self.fastq.index:
                self.fastq.at[iid, 'reverse'] = fpath
            else:
                pass

    def get_file_paths(
        self,
        result_path: str or Path,
        fasta_dir: str = None,
        fastq_dir: str = None
    ) -> None:

        if fasta_dir is None and fastq_dir is None:
            raise ValueError(
                'Please specify a directory name for FASTA or FASTQ.'
            )

        print('Getting file paths:', result_path)

        Path(result_path) if isinstance(result_path, str) else result_path

        if fasta_dir:
            self._add_fasta(result_path / fasta_dir)
        if fastq_dir:
            self._add_fastq(result_path / fastq_dir)

        self.update_iid()

    def by_iid(self) -> pandas.DataFrame:

        """ Returns a bool summary of IIDs in each process """

        exclude_attrs = ('fasta', 'fastq', 'iid')

        df = pandas.DataFrame(
            index=self.iid.iid,
            columns=[attr for attr, _ in self if attr not in exclude_attrs]
        )

        for attr, data in self:
            if attr in exclude_attrs:
                continue
            for iid in self.iid.iid:
                if iid in data.index:
                    df.at[iid, attr] = True
                else:
                    df.at[iid, attr] = False

        return df


class SurveyData(ResultData):

    def __init__(self):

        self.mlst = pandas.DataFrame()
        self.kraken = pandas.DataFrame()
        self.mash = pandas.DataFrame()
        self.kleborate = pandas.DataFrame()

        self.abricate_resistance = pandas.DataFrame()
        self.abricate_virulence = pandas.DataFrame()
        self.abricate_plasmid = pandas.DataFrame()

        self.mykrobe_phenotype = pandas.DataFrame()
        self.mykrobe_genotype = pandas.DataFrame()
        self.mykrobe_lineage = pandas.DataFrame()

    def __copy__(self):

        sd = SurveyData()
        for attr, data in self:
            setattr(sd, attr, data)
        return sd

    @property
    def empty(self):

        check = [data.empty for attr, data in self]

        if all(check):
            return True
        else:
            return False


class SketchySurvey(PoreLogger):

    """ Build sketch-associated genotype files from Pathfinder Surveys """

    def __init__(self, survey_directory: Path):

        PoreLogger.__init__(self)

        self.survey_data = SurveyData()

        self.logger.info(
            f'Parse survey directory: {survey_directory}'
        )

        self.missing = '-'

        self.survey_data.read(survey_directory)

    def construct(self, config: dict, binary: dict, merge: dict = None):

        """ Wrapper for constructing sketchy data from surveys """

        return self.construct_sketchy_data(
            config=config, binary=binary, merge=merge
        )

    def construct_sketchy_data(
        self, config: dict, binary: dict = None, merge: dict = None
    ) -> pandas.DataFrame:

        """ Construct genotype meta data for species-sketches in Sketchy

        :param config: dictionary in data: [columns] format to construct
            a dataframe with index: iids, columns: genotypes

        :param binary: dictionary of columns to treat as binary data
            where keys are columns and values are character to treat
            as missing data

        :raises: ValueError if genotype column not in selected data

        :return: dataframe of genotypes with genotype indices

        """

        genotypes = dict()
        genotype_index = None  # filled with last common index of genotypes parsed
        for attr, columns in config.items():
            data = getattr(self.survey_data, attr)
            if merge:
                for new_column, to_merge in merge[attr].items():
                    if len(to_merge) < 2:
                        raise ValueError('Merge data must contain two or more entries.')
                    try:
                        data[new_column] = data[to_merge].apply(
                            lambda x: self.missing if all(
                                # All are missing, assign missing so binary does not get confused
                                [True if y == self.missing else False for y in x]
                            ) else ';'.join(x), axis=1
                        )
                    except KeyError:
                        raise ValueError(f'Could not detect {to_merge} in data frame')

                    # Make sure name not same:
                    merge_drop = [c for c in to_merge if c != new_column]
                    data.drop(columns=merge_drop, inplace=True)

            for genotype in columns:
                if genotype not in data.columns.values:
                    raise ValueError(
                        f'Could not find {genotype} in {attr} data.'
                    )
                else:
                    genotype_data = data[genotype].tolist()
                    genotype_index = data.index
                    try:
                        _ = binary[attr]
                        try:
                            genotype_data = [
                                1 if g != self.missing else 0 for g in genotype_data
                            ]
                        except KeyError:
                            pass  # TODO: catch error?
                    except KeyError:
                        pass

                    print(genotype, len(genotype_data), len(genotype_index))

                    genotypes[genotype] = genotype_data

        if genotype_index is None:
            raise ValueError('No genotypes could be parsed.')

        return pandas.DataFrame(genotypes, genotype_index)

