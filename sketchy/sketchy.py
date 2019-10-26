"""
Rewrite of the Sketchy prediction module - cleaner, faster, more robust!
"""

import pandas
import pysam

import seaborn as sns
import matplotlib.pyplot as plt

from statistics import mode

from pathfinder.pipelines.data import SurveyData
from sketchy.utils import PoreLogger, run_cmd
from pathlib import Path

from numpy import reshape


class Sketchy(PoreLogger):

    """ Sketchy prediction API """

    def __init__(
        self,
        fastx: Path,
        reads: int = 10,
        tmp: Path = Path.cwd() / '.tmp'
    ):

        """ Sketchy prediction API

        :param fastx: Path to (compressed) input read file Fast{a,q}
        :param reads: Number of reads to extract for prediction
        :param tmp: Path to temporary directory for processing output

        """

        PoreLogger.__init__(self)

        self.fastx = fastx
        self.reads = reads

        self.tmp = tmp
        self.tmp.mkdir(parents=True, exist_ok=True)

        # Interim updates of sum of shared hashes
        self.ssh: pandas.DataFrame or None = None

    def extract_reads(self) -> Path:

        """ Extract reads for prediction with Sketchy

        Reads are extracted into temporary directory as single enumerated
        Fasta files. Creates a list of the output files for input to MASH.

        :return Path to list of output read files `mash.in` for input to MASH

        """

        self.logger.info('Extracting reads for prediction ...')

        mash_list = self.tmp / 'mash.in'
        mash_handle = mash_list.open('w')

        with pysam.FastxFile(self.fastx) as fin:
            for i, read in enumerate(fin):
                fasta = self.tmp / f'{i}.fa'
                with fasta.open('w') as fout:
                    fout.write(f">{read.name}\n{read.sequence}")
                mash_handle.write(
                    str(fasta.absolute()) + '\n'
                )
                if self.reads is not None and i == self.reads:
                    break

        return mash_list

    def compute_shared_hashes(
        self,
        reference: Path,
        query: Path,
        threads: int = 4,
        p_value: float = 1e-06
    ) -> Path:

        """ Compute shared hashes with MASH

        :param reference: Path to species-wide reference sketch created with MASH
        :param query: Path to list of extracted reads to process with MASH
        :param threads: Threads used to compute shared hashes with MASH
        :param p_value: Filter on p-value computed with MASH

        :return Path to file containing p-value filtered results from
            computing shared hashes with MASH

        """

        self.logger.info('Computing shared hashes with MASH')

        result = self.tmp / 'mash.tsv'
        run_cmd(
            f"mash dist -l -v {p_value} -p {threads}"
            f" {reference} {query} > {result}",
            shell=True
        )

        return result

    def compute_ssh(
        self,
        file: Path,
        ranks: int = 10
    ) -> pandas.DataFrame:

        """ Compute the sum of shared hashes from output of MASH

        Ranked extractions are essential at this stage to reduce the size
        of the sum of shared hashes tables at each read for large species
        sketches such as for Staphyloccocus aureus.

        :param file: Path to output file containing shared hashes from MASH

        :param ranks: Extract the highest ranking sum of shared hashes for
            evaluation on the ranked sums of sum of shared hashes

        :return Ranked sum of shared hashes against the sketch at each read

        """

        self.logger.info('Computing sum of shared hashes with Sketchy')

        converters = {
            0: lambda x: str(Path(x).name),
            1: lambda x: int(Path(x).stem),
            4: lambda x: int(x.split('/')[0])
        }

        df = pandas.read_csv(
            file,
            sep='\t',
            header=None,
            index_col=0,
            usecols=[0, 1, 4],
            names=['id', 'read', 'shared'],
            converters=converters
        )

        ssh_data = dict()
        for read, data in df.groupby('read'):

            read = int(read)
            data = data.rename(
                columns={'shared': 'ssh'}
            )
            data = data.drop('read', axis=1)

            if self.ssh is None:
                self.ssh = data
            else:
                # Sum of shared hashes, updates SSH:
                self.ssh = self.ssh.add(data, fill_value=0)

            # Sort the current SSH scores:
            self.ssh = self.ssh.sort_values('ssh', ascending=False)

            # Select the top ranking genomes by sum of shared hashes
            top_ssh = self.ssh[:ranks]

            # Add columns for evaluation plots
            top_ssh = top_ssh.assign(
                read=[read for _ in top_ssh.ssh],
                rank=[i for i in range(ranks)]
            )

            ssh_data[read] = top_ssh

        # Check if dictionary is empty, this should not be the case
        if not bool(ssh_data):
            raise ValueError(
                'No reads were classified [Sketchy Error 1]'
            )

        # First check that the first read was classified:
        if 0 not in ssh_data.keys():
            self.logger.debug(
                'First read not classified, substituting null values '
                'for initiation of sum of shared hashes summary'
            )
            # If not, take the first classified read
            # and set the top of sum of shared hashes to 0:
            ssh_data[0] = self._init_ssh(ssh_data)

        # Continuity clause, if read classification is missing due to
        # filter step on p-value or failed classification, duplicate
        # previous top sum of shared hashes to emulate null values:

        ssh = list()
        current = ssh_data[0]
        for r in range(self.reads):
            if r in ssh_data.keys():
                current = ssh_data[r]
            else:
                current = current.assign(
                    read=[r for _ in current.ssh]
                )
            ssh.append(current)

        return pandas.concat(ssh).sort_values(['read', 'rank'])

    @staticmethod
    def _init_ssh(ssh_data: dict) -> pandas.DataFrame:

        """ Helper method to initiate data entry

        Ranked null value sum of shared hashes if the first read did not contain
        any matches with the sketch; essential for continuity clause
        in computation of sum of shared hashes.

        :param ssh_data: Dictionary of read: ranked sum of shared hashes data
            frames to extract first classified entry and replace shared hashes

        :return First classified data frame with read and shared hashes replaced
            with 0 (for first read, and for null values on unclassified read)

        """

        first_classified = sorted(ssh_data.keys())[0]
        init_ssh = ssh_data[first_classified]

        init_ssh = init_ssh.assign(
            shared=[0 for _ in init_ssh.shared],
            read=[0 for _ in init_ssh.shared]
        )

        return init_ssh


class SketchyEvaluation(PoreLogger):

    def __init__(
        self,
        ssh_data: Path,
        feature_data: Path,
        top: int = 5,
        limit: int or None = 1000
    ):

        PoreLogger.__init__(self)

        self.ssh_data: pandas.DataFrame = pandas.read_csv(
            ssh_data, sep='\t', index_col=0, header=0, dtype={
                'id': str, 'ssh': float, 'read': int, 'rank': int,
            }
        )

        self.feature_data: pandas.DataFrame = pandas.read_csv(
            feature_data, sep='\t', index_col=0, header=0, dtype=str
        )

        # Make sure the index is named correctly:
        self.feature_data.index.names = ['id']

        # Merge sum of shared hashes and feature data:
        self.ssh_features = self.ssh_data \
            .join(self.feature_data, how='inner') \
            .sort_values(['read', 'rank'])

        # Settings
        self.top = top
        self.limit = limit
        self.breakpoints: dict = dict()

        self.ranks = sorted(
            self.ssh_features['rank'].unique().tolist()
        )[-1]+1

        self.reads: int = self.ssh_features.read.tolist()[-1]+1

        if self.limit is None:
            self.limit = self.reads

        # Other
        self.first_breakline: bool = False  # plot first detection breakpoint
        self.na_color: str = "#d9d9d9"      # background color heatmap

    def evaluate(
        self,
        features: list or None = None,
        stable: int = 500,
        color: str = 'YlGnBu',
        fout: Path = Path.cwd() / 'sketchy.png',
    ):

        # Something odd with colors, need reverse palettes:
        if not color.endswith('_r'):
            color += '_r'

        # Feature selection for evaluation and plotting
        if features is None:
            features: list = self.feature_data.columns.values
        else:
            # Check presence of user selected features:
            for f in features:
                if f not in self.feature_data.columns.values:
                    self.logger.info(f'Ignore feature not present in data: {f}')

        # Setup plots:
        number_features = len(features)

        fig, axes = plt.subplots(
            nrows=number_features, ncols=2, figsize=(
                14.0, number_features * 4.5
            )
        )

        if axes.ndim == 1:
            axes = reshape(
                axes, (-1, 2)
            )

        fig.subplots_adjust(hspace=0.8)

        # Feature evaluation and plot loop:
        for i, feature in enumerate(features):
            sssh, top_features = self.compute_sssh(feature=feature)

            self.breakpoints[feature] = self.define_breakpoints(
                feature=feature, top_features=top_features, stable_size=stable
            )

            self.plot_heatmap(
                feature=feature, color=color, ax=axes[i, 0]
            )

            self.plot_sssh(
                feature=feature, color=color, ax=axes[i, 1]
            )

        plt.tight_layout()
        fig.savefig(fout)

    def compute_sssh(self, feature: str) -> (pandas.DataFrame, list):

        """ Compute the sum of sums of shared hashes for a feature

        Sum of sums of shared hashes per feature, where the feature is
        the a genotype column in the species sketch associated data file.

        :param feature: column name in sketch genotype file

        :return a dataframe with sum of sums of shared hashes per read for
            the given feature and the top

        """

        self.logger.debug(
            f'Compute sum of sums of shared hashes for feature: {feature}'
        )

        top_features = []
        sssh_per_read = []
        for read, read_group in self.ssh_features.groupby('read'):
            read = int(read)

            # For each read sum the grouped feature's sums of shared hashes
            sssh = read_group.groupby(feature).sum().reset_index() \
                .sort_values('ssh', ascending=False).reset_index()

            sssh = sssh.rename(
                columns={'ssh': 'sssh'}
            )
            sssh = sssh.assign(
                rank=[i for i, _ in enumerate(sssh.sssh)],
                read=[read for _ in sssh.sssh]
            )

            top_features.append(str(
                sssh.loc[0, feature]
            ))
            sssh_per_read.append(sssh)

            if read == self.limit:
                break  # speedup

        df = pandas.concat(sssh_per_read)
        df = df.reset_index(drop=True).drop(df.columns[0], axis=1)

        return df, top_features

    def define_breakpoints(
        self, feature: str, top_features: list, stable_size: int = 500
    ) -> (int or None, int or None):

        """ Detect breakpoints of first and stable detection

        Breakpoints are computed for the most frequent feature value
        across the highest ranking sum of sums of shared hashes per read.

        :param feature: feature to evaluate breakpoint on (explicit)

        :param top_features: top ranking feature values per read
            by highest ranking sum of sums of shared hashes

        :param stable_size: number of reads for which the top
            ranking feature must be stable (continuous)

        :return: first and stable detection breakpoints

        """

        top_feature = mode(top_features)

        try:
            first = top_features.index(top_feature)+1
        except ValueError:
            first = None

        stable = None
        for i, p in enumerate(top_features):
            if p == top_feature:
                if i+stable_size > len(top_features):
                    stable = None
                else:
                    uniform = set(
                        top_features[i:i+stable_size]
                    )
                    if len(uniform) == 1 and list(uniform)[0] == top_feature:
                        stable = i+1
                        break

        self.logger.info(
            f'{feature.capitalize()} >> {top_feature} << breakpoints: '
            f'{first} (first) / {stable} (stable)'
        )

        return first, stable

    def plot_sssh(
        self,
        feature: str,
        ax: plt.axes = None,
        color: str = 'YlGnBu'
    ):

        """ Plot sum of sums of shared hashes by read for feature values

        Computation of sum of sums of shared hashes is included in the
        plotting operation using the `sum` estimator to aggregate the
        sum of sums of shared hashes for each feature value.

        :param feature:
        :param ax:
        :param color:
        :return:
        """

        palette = sns.color_palette(color, self.top)

        # Get top feature value strings by total sum of shared hashes to plot
        top_features_values = self._get_top_feature_values(feature)

        # Make sure the feature is a string for subsetting
        self.ssh_features[feature] = self.ssh_features[feature].astype(str)

        # Plot only top features
        df = self.ssh_features[
            self.ssh_features[feature].isin(top_features_values)
        ]

        # Setup plot values
        df.rename(
            columns={feature: feature.capitalize()}
        )
        df = df[df['read'] <= self.limit]
        palette = palette[:len(top_features_values)]

        # Plot breaklines and sum of sums of shard hashes plot

        self._add_breaklines(feature=feature, ax=ax)

        p2 = sns.lineplot(
            data=df, x='read', y='ssh', hue=feature, hue_order=top_features_values,
            ci=None, estimator='sum', ax=ax, palette=palette
        )

        # Legend and label settings

        legend = ax.legend()
        legend.texts[0].set_text(
            feature.capitalize()
        )

        p2.set_ylabel(
            'Sum of ranked sums of shared hashes\n', fontsize=9
        )
        p2.set_xlabel('\nReads', fontsize=9)

        p2.tick_params(labelsize=6)

        return p2

    def plot_heatmap(
        self,
        feature: str,
        top: int = 5,
        ax: plt.axes = None,
        color: str = 'YlGnBu'
    ):

        palette = sns.color_palette(color, top)

        top_features = self._get_top_feature_values(feature)

        df = self.ssh_features.copy().reset_index()

        self.logger.debug(
            f'Sum of shared hashes dataframe contains {len(df)} entries'
        )

        lineage_keys = {l: i for i, l in enumerate(top_features)}

        df = df.assign(
            vals=[lineage_keys.get(l, None) for l in df[feature]]
        )

        hm = df.pivot('rank', 'read', 'vals')
        hm = hm[hm.columns].astype(float)

        palette = palette[:len(top_features)]

        p1 = sns.heatmap(
            hm.iloc[:self.ranks, :self.limit], linewidths=0, cbar=False,
            ax=ax, cmap=palette
        )

        p1.set_facecolor(self.na_color)

        xticks, yticks = self._get_ticks()

        p1.set_xticks(xticks)
        p1.set_xticklabels(xticks, rotation='vertical')

        p1.set_xlabel('\nReads', fontsize=9)
        p1.set_ylabel('Ranked sum of shared hashes\n', fontsize=9)

        p1.set_yticks(yticks)
        p1.set_yticklabels(yticks)

        p1.tick_params(axis='both', which='major', labelsize=6)
        p1.tick_params(length=1, width=0.5)

        return p1

    def _get_ticks(self):

        if self.reads <= 200:
            xticks = [i for i in range(0, self.reads + 1, 20)]
        elif 200 < self.reads <= 500:
            xticks = [i for i in range(0, self.reads + 1, 50)]
        elif 500 < self.reads <= 1500:
            xticks = [i for i in range(0, self.reads + 1, 200)]
        elif 1500 < self.reads <= 5000:
            xticks = [i for i in range(0, self.reads + 1, 1000)]
        elif 5000 < self.reads <= 15000:
            xticks = [i for i in range(0, self.reads + 1, 3000)]
        else:
            xticks = [i for i in range(0, self.reads + 1, 5000)]

        if self.ranks > 10:
            yticks = [i for i in range(0, self.ranks + 1, 10)]
        else:
            yticks = [i for i in range(0, self.ranks + 1)]

        return xticks, yticks

    def _add_breaklines(self, feature: str, ax) -> None:

        try:
            bpoints = self.breakpoints[feature]
            if bpoints[0] is not None and self.first_breakline:
                ax.axvline(
                    x=bpoints[0], linewidth=1, linestyle='--', color='black'
                )

            if bpoints[1] is not None:
                ax.axvline(
                    x=bpoints[1], linewidth=1, linestyle='-', color='black'
                )
        except KeyError:
            self.logger.debug(f'Breakpoints for {feature} do not exist.')

    def _get_top_feature_values(self, feature: str) -> [str]:

        return self.ssh_features.groupby(by=feature).sum().ssh \
                 .sort_values(ascending=False)[:self.top] \
                 .index.tolist()  # TODO: check to remove NA

    @staticmethod
    def compute_preference_scores(
        sssh: pandas.DataFrame, feature: str
    ) -> None:

        """ Compute a preference score of sums of sums of shared hashes

        Analogous to Brinda et al. (2019) lineage score:

            LS = 2f/(f+t)-1

        where f and t denote weights of best matches, here these are the
        two feature values with the highest ranking total sums of
        sums of shared hashes or sum of shared hashes
        across the prediction data.

        :param sssh: dataframe for sssh preference scores
        :param feature: feature to compute the preference score for

        """

        feature = sssh.columns.values[0]

        pass


class SketchySurvey(PoreLogger):

    """ Build sketch-associated genotype files from Pathfinder Surveys """

    def __init__(self, survey_directory: Path):

        PoreLogger.__init__(self)

        self.survey = SurveyData()

        self.logger.info(
            f'Parse survey directory: {survey_directory}'
        )
        self.survey.read(survey_directory)

    def construct(self, config: dict, binary: dict):

        return self.survey.construct_sketchy_data(
            config=config, binary=binary
        )






