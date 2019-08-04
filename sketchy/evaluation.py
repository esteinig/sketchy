import pandas
import numpy as np
import scipy.stats
import seaborn as sns

from tqdm import tqdm
from sketchy.minhash import MashScore
from pathlib import Path
from cycler import cycler
from matplotlib import pyplot as plt
from sketchy.utils import PoreLogger

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
        true_lineage: str or None = None,
        true_genotype: str or None = None,
        true_resistance: str or None = None,
        palette: str = None,
        primary_color: str = "#41b6c4",
        secondary_color: str = "#7fcdbb",
        sequential: bool = True,
        sketch_data: pandas.DataFrame or Path = None
    ):

        self.indir = indir  # Temporary output directory of predict + --keep

        self.limit = limit
        self.outdir = outdir

        self.outdir.mkdir(parents=True, exist_ok=True)

        self.true_lineage = true_lineage
        self.true_susceptibility = true_resistance
        self.true_genotype = true_genotype

        self.false_color = "#d9d9d9"

        if palette in ('blue', 'red', 'green', 'orange'):
            self.true_color = BREWER[palette][-3]
            self.lineage_color = BREWER[palette][3]
        else:
            self.true_color = primary_color
            self.lineage_color = secondary_color

        self.logger = PoreLogger().logger

        self.top = top

        if isinstance(sketch_data, Path):
            self.sketch_data = pandas.read_csv(
                sketch_data, sep='\t', index_col=0
            )
        else:
            self.sketch_data = sketch_data

        self.reads: int = 0

        self.top_ssh: pandas.DataFrame = self._parse_hashes(
            sequential=sequential
        )

        self.logger.info(
            f'{len(self.top_ssh.index.unique())} unique genome hits were '
            f'contained in the top {self.top} ranked sum of shared hashes'
            f'across {len(self.top_ssh.read.unique())} reads'
        )
        self.logger.info(
            f'Genomes span {len(self.top_ssh.lineage.unique())} lineages, '
            f'{len(self.top_ssh.susceptibility.unique())} resistance profiles, '
            f'and {len(self.top_ssh.genotype.unique())} genotypes'
        )

        validate = [
            self.true_lineage, self.true_susceptibility, self.true_genotype
        ]

        if validate.count(None) == len(validate):
            self.logger.info(
                f'Generate lineage plots for evaluation on a total of '
                f'{len(self.top_ssh.read.unique())} reads'
            )
        else:

            self.logger.info(
                f'Validation: {self.true_lineage} :: '
                f'{self.true_susceptibility} :: '
                f'{self.true_genotype}'
            )
            self.top_ssh = self._assign_truth(self.top_ssh)

            # Parse complete data and assign truth for validation,
            # this may take some time
            self.top_ssh_all: pandas.DataFrame = self._parse_all_hashes()

            self.logger.info(
                'Assigning truth data, this may take some time...'
            )
            self.top_ssh_all = self._assign_truth(
                self.top_ssh_all, category=True
            )

        self.breakpoint_detection = 0
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

    def create_validation_hitmap(self, ranks: int = 50, ax=None):

        self.logger.info(
            f'Generate validation hitmap: {self.outdir / "heatmap.pdf"}'
        )

        self.top_ssh_all.reset_index(inplace=True)

        hm = self.top_ssh_all.pivot('rank', 'read', 'truth')

        hm = hm[hm.columns].astype(float)

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

        self.logger.info(
            f'Generate lineage race plot: {self.outdir / "race_plot.pdf"}'
        )

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

        self.logger.info(
            f'Generate trait concordance plot: {self.outdir / "race_plot_mean_95.pdf"}'
        )

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

        p2.set_ylabel('Mean sum of shared hashes', fontsize=10)
        p2.set_xlabel('Read', fontsize=10)

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

    def create_lineage_hitmap(self, top: int = 10, ax=None):

        self.logger.debug(
            f'Generate lineage hitmap'
        )

        self.top_ssh.lineage = self.top_ssh.lineage.astype(str)
        self.top_ssh.shared = self.top_ssh.shared.astype(float)

        top_lineages = self.top_ssh.groupby(by='lineage') \
            .sum().shared.sort_values(ascending=False)[:top]

        top_lineages = top_lineages.index.tolist()
        top_lineages = [t for t in top_lineages if t != 'nan']

        df = self.top_ssh.copy().reset_index()

        self.logger.debug(
            f'Sum of shared hashes dataframe contains {len(df)} entries'
        )

        lineage_keys = {l: i for i, l in enumerate(top_lineages)}

        df = df.assign(
            vals=[lineage_keys.get(l, None) for l in df.lineage]
        )

        hm = df.pivot('rank', 'read', 'vals')

        hm = hm[hm.columns].astype(float)

        p1 = sns.heatmap(
            hm.iloc[:self.top, :], linewidths=0, cbar=False, ax=ax,
            cmap=sns.hls_palette(top)
        )

        p1.set_facecolor(self.false_color)

        if self.reads <= 200:
            xticks = [i for i in range(0, self.reads+1, 20)]
        elif 200 < self.reads <= 500:
            xticks = [i for i in range(0, self.reads+1, 50)]
        elif 500 < self.reads <= 1500:
            xticks = [i for i in range(0, self.reads+1, 500)]
        elif 1500 < self.reads <= 5000:
            xticks = [i for i in range(0, self.reads+1, 1000)]
        elif 5000 < self.reads <= 15000:
            xticks = [i for i in range(0, self.reads + 1, 3000)]
        else:
            xticks = [i for i in range(0, self.reads+1, 5000)]

        p1.set_xticks(xticks)
        p1.set_xticklabels(xticks, rotation='vertical')

        p1.set_xlabel('\nReads', fontsize=10)
        p1.set_ylabel('Ranked sum of shared hashes\n', fontsize=10)

        if self.top > 10:
            yticks = [i for i in range(0, self.top+1, 10)]
        else:
            yticks = [i for i in range(0, self.top+1)]

        p1.set_yticks(yticks)
        p1.set_yticklabels(yticks)

        p1.tick_params(axis='both', which='major', labelsize=6)
        p1.tick_params(length=1, width=0.5)

        if self.breakpoint_stable:
            plt.axvline(x=self.breakpoint_stable, linewidth=1, color='black')

        plt.close()

        return p1

    def create_lineage_plot(
            self, top: int = 10, ax=None
    ):

        self.logger.debug(
            f'Generate lineage total sum plot'
        )

        # Changed from mean to sum() here:
        top_lineages = self.top_ssh.groupby(by='lineage')\
            .sum().shared.sort_values(ascending=False)[:top].index.tolist()

        df = self.top_ssh[self.top_ssh.lineage.isin(top_lineages)]

        df.rename(columns={'lineage': 'Lineage'})

        p3 = sns.lineplot(
            data=df, x='read', y='shared', hue='lineage', hue_order=top_lineages,
            ci=None, estimator='sum', ax=ax,
            palette=sns.hls_palette(top)
        )

        legend = ax.legend()
        legend.texts[0].set_text("Lineage")

        p3.set_ylabel(
            'Sum total of ranked sums of shared hashes\n', fontsize=10
        )
        p3.set_xlabel('\nRead', fontsize=10)

        p3.tick_params(labelsize=6)

        if self.breakpoint_detection:
            plt.axvline(
                x=self.breakpoint_detection, linewidth=1,
                linestyle='--', color='black'
            )

        if self.breakpoint_stable:
            plt.axvline(x=self.breakpoint_stable, linewidth=1,
                        color='black')

        plt.close()

        return p3

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

        ms = MashScore()

        if sequential:

            self.logger.info(
                f'Parse read hashes output in {self.indir}'
            )

            hash_dfs = []
            for i, fpath in enumerate(sorted(
                self.indir.glob('*.counts.*'),
                key=lambda x: int(
                    x.name.split('.')[-1]
                )
            )[:self.limit]):

                df = pandas.read_csv(
                    fpath, sep="\t", index_col=0, nrows=self.top
                )  # [:self.top]

                n = int(fpath.name.split('.')[-1])

                self.logger.debug(f'Process Sketchy score at read: {n}')

                df['read'] = [n for _ in df.shared]

                hash_dfs.append(df)
                self.reads = i+1

            dd = pandas.concat(hash_dfs)  #, sort=False)

        else:

            reads = sorted([
                int(
                    fpath.name.split('.')[-1]
                ) for fpath in self.indir.glob('read.*')
            ])

            read_range = list(
                range(1, reads[-1]+1)
            )
            self.logger.debug(
                f'Detected last computed read @ {reads[-1]};'
                f' set read range: 1 - {reads[-1]}'
            )

            # Check if all reads present:
            nohit = sum(1 for i in read_range if i not in reads)

            self.logger.info(
                f'Detected {nohit} reads that did not contain hits'
                f' against the sketch'
            )

            self.logger.info(
                f'Compute sum of shared hashes ...'
            )

            # Seed complete index as baseline
            hash_top = pandas.DataFrame(
                index=self.sketch_data.index.tolist(),
                data=dict(
                    shared=[0 for _ in self.sketch_data.index]
                )
            )

            hash_dfs = []
            for i in tqdm(
                read_range[:self.limit]
            ):
                self.logger.debug(
                    f'Computing sum of shared hashes @ read {i}'
                )
                read_output = self.indir / f'read.{i}'

                if read_output.exists():
                    df = pandas.read_csv(
                        read_output, sep='\t'
                    ).set_index('uuid')

                    # Make sure complete index of sketch is
                    # used in SSH computation:
                    if ms.inter.empty:
                        ms.inter = hash_top

                    ms.compute_ssh(
                        mash=df,
                        sequential=True,
                        printout=False
                    )

                    hash_top = ms.inter[:self.top]
                    hash_top = hash_top.assign(
                        read=[i for _ in hash_top.shared],
                        rank=[i for i in range(self.top)]
                    )

                    hash_dfs.append(hash_top)
                else:
                    # Continuity clause:

                    # If read is entirely missing, update read of previous:
                    hash_top = hash_top.assign(
                        read=[i for _ in hash_top.shared],
                        rank=[i for i, _ in enumerate(hash_top.shared)]
                    )

                    hash_dfs.append(
                        hash_top[:self.top]  # Make sure only top ranking
                    )  # If first read missing (no hits) take top random genomes

                self.reads = i+1

            dd = pandas.concat(hash_dfs, sort=False)
            dd = dd.join(self.sketch_data, how='inner')  # Orders by index

        if self.limit is None:
            self.limit = self.reads

        return dd

    def _parse_all_hashes(self):
        """ Parse the hashes of the unique genomes in the top reads
        This is required because these may drop out and not be represented
        in the top matches at higher read numbers.
        """

        self.logger.info(
            f'Parse unique genome hash match outputs in {self.indir} ...'
        )

        hash_dfs = []
        for i, fpath in enumerate(sorted(
                self.indir.glob('*.counts.*'),
                key=lambda x: int(
                    x.name.split('.')[-1]
                )
        )):

            if self.limit is not None and i >= self.limit:
                break

            self.logger.debug(f'Process Sketchy score at read: {i}')

            top_hits = self.top_ssh.index.unique()

            df = pandas.read_csv(
                fpath, sep="\t", index_col=0,
            )

            df = df[
                df.index.isin(top_hits)
            ]

            df['read'] = [i for _ in df.shared]
            df['rank'] = [i for i in range(len(df))]

            hash_dfs.append(df)

        return pandas.concat(hash_dfs)









