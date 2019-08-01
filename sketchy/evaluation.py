import pandas
import numpy as np
import scipy.stats
import seaborn as sns

from tqdm import tqdm
from pathlib import Path
from cycler import cycler
from matplotlib import pyplot as plt

BREWER = {
  'blue': [
      "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0",
      "#225ea8", "#253494", "#081d58"
  ],
  'green': [
      "#f7fcb9", "#d9f0a3",	"#addd8e",	"#78c679", "#41ab5d",
      "#238443", "#006837",	"#004529"
  ],
  'red': [
      "#e7e1ef", "#d4b9da", "#c994c7",	"#df65b0", "#e7298a",
      "#ce1256", "#980043", "#67001f"
  ],
  'orange': [
      "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a",
      "#e31a1c", "#bd0026",	"#800026"
  ]
}


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array([v for v in data if v is not None])
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


class SampleEvaluator:
    """ Base access to evaluation and plotting for pre-print analysis """

    def __init__(
        self,
        indir: Path,
        outdir: Path = None,
        limit: int = 1000,
        top: int = 10,
        true_lineage: str = '8',
        true_genotype: str = 'nan-nan-blaZ-nan-nan-nan-nan-nan-nan-tetL-nan-nan',
        true_resistance: str = 'RSSSSSSSRRRS',
        palette: str = None,
        primary_color: str = "#41b6c4",
        secondary_color: str = "#7fcdbb"
    ):

        self.indir = indir  # Temporary output directory of predict + --keep

        self.limit = limit
        self.outdir = outdir

        self.outdir.mkdir(parents=True, exist_ok=True)

        self.true_lineage = true_lineage  # ACTT 243, Zymo 9
        self.true_susceptibility = true_resistance  # ACTT SSSSSSSSSSSS Zymo SRSSSSSSRSSS
        self.true_genotype = true_genotype # ACTT 'nan-nan-nan-nan-nan-nan-nan-nan-nan-nan-nan-nan'

        self.false_color = "#d9d9d9"

        if palette in ('blue', 'red', 'green', 'orange'):
            self.true_color = BREWER[palette][-3]
            self.lineage_color = BREWER[palette][3]
        else:
            self.true_color = primary_color
            self.lineage_color = secondary_color

        print(
            f'Truth: {self.true_lineage} :: '
            f'{self.true_susceptibility} :: '
            f'{self.true_genotype}'
        )

        self.top: int = top
        self.reads: int = 0

        self.top_ssh: pandas.DataFrame = self._parse_hashes()
        self.top_ssh_all: pandas.DataFrame = self._parse_all_hashes()

        print(
            f'There are {len(self.top_ssh.index.unique())} unique genomes '
            f'hit in the top {self.top} of {len(self.top_ssh.read.unique())} '
            f'reads including {len(self.top_ssh.lineage.unique())} lineages, '
            f'{len(self.top_ssh.susceptibility.unique())} resistance profiles, '
            f'and {len(self.top_ssh.genotype.unique())} genotypes.'
        )

        self.top_ssh = self._assign_truth(self.top_ssh)

        match_count = self.top_ssh.groupby(
            self.top_ssh.index
        ).apply(self._sum_matches).sum()

        print(f'There were {match_count} unique matches on lineage and traits in'
              f' the top {self.top} ssh-matches over {self.reads} reads.')

        self.top_ssh_all = self._assign_truth(
            self.top_ssh_all, category=True
        )
        self.breakpoint_detection = 0  #  self.to[df['A'] == 5].index.item()
        self.breakpoint_stable = 0

    @staticmethod
    def _sum_matches(data):
        """ Joint lineage / genotype / susceptibility match """
        d = data['truth'].values[0]
        if d == 2:
            return 1
        else:
            return 0

    def _assign_truth(self, df, name: bool = False, category: bool = True):

        grouped = df.groupby(by=df.index)

        df['truth'] = [
            'None' for _ in df.shared
        ]

        for k, v in grouped:
            _df = grouped.get_group(k)
            df.at[
                _df.index[0], 'truth'
            ] = self._get_truth_color(_df, category=category, name=name)

        return df

    def create_timeline_hitmap(self, ranks: int = 50, ax=None):

        self.top_ssh_all.reset_index(inplace=True)

        hm = self.top_ssh_all.pivot('rank', 'read', 'truth')

        hm = hm[hm.columns].astype(float)
        hm2 = hm.iloc[:ranks, :]

        p1 = sns.heatmap(
            hm.iloc[:ranks, :], linewidths=0, cbar=False, ax=ax,
            cmap=[self.false_color, self.lineage_color, self.true_color]
        )

        if self.reads <= 200:
            xticks = [i for i in range(0, self.reads+1, 20)]
        elif 200 < self.reads <= 500:
            xticks = [i for i in range(0, self.reads+1, 50)]
        elif 500 < self.reads <= 1500:
            xticks = [i for i in range(0, self.reads+1, 500)]
        else:
            xticks = [i for i in range(0, self.reads+1, 5000)]

        p1.set_xticks(xticks)
        p1.set_xticklabels(xticks, rotation='horizontal')

        p1.set_xlabel('Reads', fontsize=8)
        p1.set_ylabel('', fontsize=8)

        if ranks > 10:
            yticks = [i for i in range(0, ranks+1, 10)]
        else:
            yticks = [i for i in range(0, ranks+1)]

        print(yticks)

        # Rank based index from 1
        yticks[0] = 1

        p1.set_yticks(yticks)
        p1.set_yticklabels(yticks)

        p1.tick_params(axis='both', which='major', labelsize=6)
        p1.tick_params(length=1, width=0.5)

        if self.breakpoint_stable:
            plt.axvline(x=self.breakpoint_stable, linewidth=1, color='black')

        p1.get_figure().savefig(
            f'{self.outdir / "heatmap.pdf"}',
            figsize=(11.0, 7.0)
        )
        plt.close()

        return p1

    def create_race_plot(self, ax=None):

        # Re-assign truth as colors for plotting
        self.top_ssh_all = self._assign_truth(
            self.top_ssh_all, name=False, category=False
        )

        colors = [
            self.top_ssh_all.at[idx, 'truth']
            for idx in self.top_ssh_all.index.unique().values
        ]

        # Re-assign truth as names for plotting
        self.top_ssh_all = self._assign_truth(
            self.top_ssh_all, name=True, category=False
        )

        df = self.top_ssh_all.rename(
            {'truth': 'Concordance'}, axis=1
        )

        plt.rc(
            'axes', prop_cycle=(cycler('color', colors))
        )

        p1 = sns.lineplot(
            data=df, x='read', y='shared', hue='uuid', style='Concordance',
            legend=False, estimator=None, ci=None, lw=0.8, ax=ax
        )

        p1.set_ylabel('Sum of shared hashes', fontsize=8)
        p1.set_xlabel('Read', fontsize=8)

        p1.tick_params(labelsize=6)

        if self.breakpoint_detection:
            plt.axvline(
                x=self.breakpoint_detection, linewidth=1,
                linestyle='--', color='black'

            )

        if self.breakpoint_stable:
            plt.axvline(x=self.breakpoint_stable, linewidth=1, color='black')

        p1.get_figure().savefig(
            f'{self.outdir / "race_plot.pdf"}',
            figsize=(11.0, 7.0)
        )
        plt.close()

        return p1

    def create_concordance_plot(self, ax=None):

        df = self.top_ssh_all.rename(
            {'truth': 'Concordance'}, axis=1
        )

        p2 = sns.lineplot(
            data=df, x='read', y='shared', hue='Concordance',
            ci=None, estimator='mean', ax=ax,
            palette=sns.color_palette(
                [self.true_color, self.false_color, self.lineage_color], len(
                    df.Concordance.unique()
                )
        ))

        p2.set_ylabel('Mean sum of shared hashes', fontsize=8)
        p2.set_xlabel('Read', fontsize=8)

        p2.tick_params(labelsize=6)

        if self.breakpoint_detection:
            plt.axvline(
                x=self.breakpoint_detection, linewidth=1,
                linestyle='--', color='black'

            )

        if self.breakpoint_stable:
            plt.axvline(x=self.breakpoint_stable, linewidth=1, color='black')

        p2.get_figure().savefig(
            f'{self.outdir / "race_plot_mean_95.pdf"}',
            figsize=(11.0, 7.0)
        )
        plt.close()

        return p2

    def create_lineage_plot(self, ax=None):


        # Select 10 lineages ranked by mean SSH

        top_lineages = self.top_ssh_all.groupby(by='lineage')\
            .mean('shared').sort_values('shared', ascending=False)

        print('Test', top_lineages)

        col = sns.color_palette("muted")


        # p3 = sns.lineplot(
        #     data=df, x='read', y='shared', hue='lineage',
        #     ci=None, estimator='mean', ax=ax,
        #     palette=sns.color_palette(
        #         len(
        #             df.Concordance.unique()
        #         )
        #     ))
        #
        # p3.set_ylabel('Mean sum of shared hashes', fontsize=8)
        # p3.set_xlabel('Read', fontsize=8)
        #
        # p2.tick_params(labelsize=6)
        #
        # if self.breakpoint_detection:
        #     plt.axvline(
        #         x=self.breakpoint_detection, linewidth=1,
        #         linestyle='--', color='black'
        #
        #     )
        #
        # if self.breakpoint_stable:
        #     plt.axvline(x=self.breakpoint_stable, linewidth=1,
        #                 color='black')
        #
        # p2.get_figure().savefig(
        #     f'{self.outdir / "race_plot_mean_95.pdf"}',
        #     figsize=(11.0, 7.0)
        # )
        # plt.close()
        #
        # return p2

    def _get_truth_color(
        self,
        df: pandas.DataFrame,
        name: bool = False,
        category: bool = False
    ):
        """
        Get truthy color for each attribute in the uuid-grouped dataframe.

        Used as part of loop iterating over index-grouped DataFrame
        """

        lineage = str(df.lineage.iloc[0])
        susceptibility = str(df.susceptibility.iloc[0])
        genotype = str(df.genotype.iloc[0])

        if lineage == self.true_lineage:
            color = self.true_color
            if name:
                color = 'True'
            if category:
                color = 2
            if genotype != self.true_genotype:
                color = self.lineage_color
                if name:
                    color = 'Lineage'
                if category:
                    color = 1
            if susceptibility != self.true_susceptibility:
                color = self.lineage_color
                if name:
                    color = 'Lineage'
                if category:
                    color = 1
        else:
            color = self.false_color
            if name:
                color = 'False'
            if category:
                color = 0

        return color

    def _parse_hashes(self, sequential: bool = True) -> pandas.DataFrame:

        print(f'Parsing read hashes output in {self.indir}')
        hash_dfs = []
        if sequential:
            for i, fpath in enumerate(sorted(
                self.indir.glob('*.counts.*'),
                key=lambda x: int(
                    x.name.split('.')[-1]
                )
            )):

                df = pandas.read_csv(fpath, sep="\t", index_col=0)[:self.top]
                n = int(fpath.name.split('.')[-1])
                df['read'] = [n for _ in df.shared]
                if self.limit is not None and i >= self.limit:
                    break

                hash_dfs.append(df)
                self.reads = i+1
        else:

            inter = pandas.DataFrame()

            reads = sorted([
                int(
                    fpath.name.split('.')[-1]
                ) for fpath in self.indir.glob('read.*')
            ])

            # Check if all reads present:

            for i in range(reads[-1]):
                if i + 1 not in reads:
                    raise ValueError(
                        f'Read output directory is missing hashes at read {i}'
                    )

            for i in tqdm(reads, desc='Process read'):

                df = pandas.read_csv(
                    self.indir / f'read.{i}', sep='\t'
                ).sort_values(by='shared', ascending=False)

                if inter.empty:
                    inter = df
                else:
                    # Sum of shared hashes
                    inter = inter.add(df)

                d = inter[:self.top]
                d['read'] = [i for _ in d.shared]
                if self.limit is not None and i >= self.limit:
                    break

                hash_dfs.append(d)
                self.reads = i + 1

        if self.limit is None:
            self.limit = self.reads

        return pandas.concat(hash_dfs).sort_values(by='read')

    def _parse_all_hashes(self):
        """ Parse the hashes of the unique genomes in the top reads
        This is required because these may drop out and not be represented
        in the top matches at higher read numbers.
        """

        print(f'Parsing unique genome hash output in {self.indir}')
        hash_dfs = []
        for i, fpath in enumerate(sorted(
                self.indir.glob('*.counts.*'),
                key=lambda x: int(
                    x.name.split('.')[-1]
                )
        )):

            if self.limit is not None and i >= self.limit:
                break

            df = pandas.read_csv(fpath, sep="\t", index_col=0)
            df = df[
                df.index.isin(self.top_ssh.index.unique())
            ]
            df['read'] = [i for _ in df.shared]
            df['rank'] = [i for i in range(len(df))]

            hash_dfs.append(df)

        return pandas.concat(hash_dfs)









