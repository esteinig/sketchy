import pandas
import seaborn as sns
import itertools

from tqdm import tqdm
from sketchy.minhash import MashScore
from pathlib import Path
from matplotlib import pyplot as plt
from sketchy.utils import PoreLogger


class Sample:

    """ Base access to evaluation and plotting for pre-print analysis """

    def __init__(
        self,
        indir: Path = None,
        outdir: Path = None,
        limit: int = 1000,
        top: int = 10,
    ):

        self.top = top
        self.indir = indir
        self.limit = limit

        self.outdir = outdir
        if outdir:
            self.outdir.mkdir(parents=True, exist_ok=True)

        self.false_color = "#d9d9d9"

        self.logger = PoreLogger().logger

        self.reads: int = 0

        self.breakpoint_detection: int or None = None
        self.breakpoint_stable: int or None = None

    def _parse_hashes(self) -> pandas.DataFrame:

        """
        Parse output of sketchy predict temporary directory
         - overwritten in subclasses
        """

        pass



class SampleEvaluator(Sample):

    def __init__(
        self,
        indir: Path = None,
        outdir: Path = None,
        limit: int = 1000,
        top: int = 10,
        plot_data: pandas.DataFrame = None,
        sketch_data: pandas.DataFrame or Path = None
    ):

        super().__init__(
            indir=indir, outdir=outdir, limit=limit, top=top
        )

        if isinstance(sketch_data, Path):
            self.sketch_data = pandas.read_csv(
                sketch_data, sep='\t', index_col=0
            )
        else:
            self.sketch_data = sketch_data

        if indir is not None:
            self.top_ssh: pandas.DataFrame = self._parse_hashes()
            self._log_info()

        if plot_data is not None:
            self.top_ssh = plot_data

        if self.reads == 0:
            self.reads = self.limit

    def _log_info(self):

        self.logger.info(
            f'{len(self.top_ssh.index.unique())} unique genome hits were '
            f'contained in the top ranking sum of shared hashes ({self.top}).'
        )
        self.logger.info(
            f'Genomes span {len(self.top_ssh.lineage.unique())} lineages, '
            f'{len(self.top_ssh.susceptibility.unique())} resistance profiles, '
            f'and {len(self.top_ssh.genotype.unique())} genotypes'
        )

        self.logger.info(
            f'Generate lineage plots for evaluation on a total of '
            f'{len(self.top_ssh.read.unique())} reads'
        )

    def _parse_hashes(self) -> pandas.DataFrame:

        ms = MashScore()

        reads = sorted([
            int(
                fpath.name.split('.')[-1]
            ) for fpath in self.indir.glob('read.*')
        ])

        read_range = list(
            range(1, reads[-1] + 1)
        )
        self.logger.info(
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

            self.reads = i + 1

        dd = pandas.concat(hash_dfs, sort=False)
        dd = dd.join(self.sketch_data, how='inner')  # Orders by index

        return dd

    def sum_of_ssh(self, data):

        self.top_ssh.lineage = self.top_ssh[data].astype(str)
        self.top_ssh.shared = self.top_ssh.shared.astype(float)

        # Compute the sum of sum of shared hashes

        self.logger.info(
            f'Computing sum of sums of shared hashes'
            f' for {data.capitalize()} ...'
        )
        top_predictions = []
        i = 0
        for _, read_group in self.top_ssh.groupby('read'):
            i += 1
            data_sums = read_group.groupby(data).sum().reset_index() \
                .sort_values('shared', ascending=False).reset_index()
            top_predictions.append(
                str(data_sums.loc[0, data])
            )
        self.logger.info(f'This many reads grouped: {i}')
        return top_predictions

    def find_breakpoints(
            self, top: int = 5, data: str = 'lineage', block_size: int = 500
    ):

        top_predictions = self.sum_of_ssh(data)

        top_data = self.top_ssh.groupby(by=data) \
            .sum().shared.sort_values(ascending=False)[:top]

        topaz = str(top_data.index.tolist()[0])

        self.breakpoint_detection = top_predictions.index(topaz)+1
        self.logger.debug(
            f'Found detection breakpoint @'
            f' {self.breakpoint_detection} reads.'
        )
        self.logger.info('Computing stable detection breakpoint ...')
        for i, p in enumerate(top_predictions):
            if p == topaz:
                # Check for limit ahead:
                uniform = set(
                    top_predictions[i:i+block_size]
                )

                if len(uniform) == 1 and list(uniform)[0] == topaz:
                    self.breakpoint_stable = i
                    break

        # TODO might change to most common lineage at last evaluated read
        self.logger.debug(
            f'Found stable breakpoint @'
            f' {self.breakpoint_stable} reads.'
        )
        self.logger.info(
            f'Breakpoints for {data.capitalize()} {topaz} '
            f'@ {self.breakpoint_detection} reads (first) '
            f'and {self.breakpoint_stable} reads (stable)'
        )

    def create_hitmap(
            self, top: int = 5, data: str = 'lineage', ax=None, color='hls'
    ):

        self.logger.debug(
            f'Generate lineage hitmap [ {top}, {data}]'
        )

        self.top_ssh[data] = self.top_ssh[data].astype(str)
        self.top_ssh.shared = self.top_ssh.shared.astype(float)

        palette = sns.color_palette(color, top)

        top_lineages = self.top_ssh.groupby(by=data) \
            .sum().shared.sort_values(ascending=False)[:top]

        top_lineages = top_lineages.index.tolist()
        top_lineages = [t for t in top_lineages if t != 'nan']

        df = self.top_ssh.copy().reset_index()

        self.logger.debug(
            f'Sum of shared hashes dataframe contains {len(df)} entries'
        )

        lineage_keys = {l: i for i, l in enumerate(top_lineages)}

        df = df.assign(
            vals=[lineage_keys.get(l, None) for l in df[data]]
        )

        hm = df.pivot('rank', 'read', 'vals')

        hm = hm[hm.columns].astype(float)

        p1 = sns.heatmap(
            hm.iloc[:self.top, :self.limit], linewidths=0, cbar=False,
            ax=ax, cmap=palette
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

        p1.set_xlabel('\nReads', fontsize=9)
        p1.set_ylabel('Ranked sum of shared hashes\n', fontsize=9)

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

        return p1

    def create_lineplot(
            self, top: int = 10, data: str = 'lineage', ax=None, color='hls'
    ):

        self.logger.debug(
            f'Generate lineage total sum plot'
        )

        palette = sns.color_palette(color, top)

        # Changed from mean to sum() here:
        top_lineages = self.top_ssh.groupby(by=data)\
            .sum().shared.sort_values(ascending=False)[:top].index.tolist()

        df = self.top_ssh[self.top_ssh[data].isin(top_lineages)]

        df.rename(columns={data: data.capitalize()})

        df = df[df['read'] <= self.limit]

        p3 = sns.lineplot(
            data=df, x='read', y='shared', hue=data, hue_order=top_lineages,
            ci=None, estimator='sum', ax=ax, palette=palette
        )

        legend = ax.legend()
        legend.texts[0].set_text(data.capitalize())

        p3.set_ylabel(
            'Sum of ranked sums of shared hashes\n', fontsize=9
        )
        p3.set_xlabel('\nReads', fontsize=9)

        p3.tick_params(labelsize=6)

        if self.breakpoint_detection:
            plt.axvline(
                x=self.breakpoint_detection, linewidth=1,
                linestyle='--', color='black'
            )

        if self.breakpoint_stable:
            plt.axvline(x=self.breakpoint_stable, linewidth=1,
                        color='black')

        return p3
