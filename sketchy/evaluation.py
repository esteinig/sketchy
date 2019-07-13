from sketchy.sketchy import Sketchy
from pathlib import Path

import shutil
import pandas
import random
import sys

import numpy as np
import scipy.stats
import seaborn as sns

from tqdm import tqdm
from matplotlib import pyplot as plt

from sketchy.utils import run_cmd


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array([v for v in data if v is not None])
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


class Evaluator:
    """ Base access class of algorithm evaluation for the manuscript """

    def __init__(self, outdir: Path = None):

        self.sketchy = Sketchy()

        self.outdir = outdir

        if outdir:
            self.outdir.mkdir(parents=True, exist_ok=True)

        self.boot_prefix = 'boot'

    def __iter__(self):

        dirs = [
            path for path in self.outdir.glob('*')
            if path.is_dir() and not path.name.endswith('_tmp')
        ]

        for directory in dirs:
            yield directory, list(directory.glob(f"{self.boot_prefix}*"))

    def _identify_stages(self):
        """ First, first correct stable, isolate stable """
        pass

    def bootstrap(
        self,
        fastq: Path,
        nbootstrap: int = 10,
        sample_reads: int or None = None,
        sample_read_proportion: float or None = 1.0,
        shuffle=True
    ) -> [Path]:

        print(f'Bootstrapping reads in file: {fastq}')

        bsdir = fastq.stem
        (self.outdir / bsdir).mkdir(parents=True, exist_ok=False)

        if shuffle:
            # Required because fastq-sample otherwise respects
            # order of reads in FASTQ file - this would therefore sample reads
            # in relative time intervals cf. to the original run, which
            # defeats the point of random samples with replacement in
            # bootstraps, so reads have to be shuffled first:
            seed_sort = random.randint(0, sys.maxsize)
            randomized = Path(
                str(self.outdir / bsdir / fastq.stem) + '_random.fq'
            )
            run_cmd(
                f'fastq-sort --random --seed={seed_sort} {fastq} > {randomized}',
                shell=True
            )
            fastq = randomized

        reads = f'-n {sample_reads}' if sample_reads else f'-p {sample_read_proportion}'

        for i in tqdm(range(nbootstrap), desc='Sample bootstrap replicates:'):
            seed_bs = random.randint(0, sys.maxsize)
            boot_file = self.outdir / bsdir / f"{self.boot_prefix}{i}"
            run_cmd(
                f'fastq-sample -s {seed_bs} -r {reads}' 
                f' -o {boot_file} {fastq}',
                shell=True
            )

        return list(
            (self.outdir / bsdir).glob('*.fastq')
        )

    def predict_bootstraps(
        self,
        bsfiles,
        sketch,
        data,
        cores=4,
        reads=100
    ) -> pandas.DataFrame:

        print(f'Predicting scores on bootstrap files.')

        score_data = []
        for fq in tqdm(bsfiles, desc='Bootstrap replicates'):
            replicate_name = fq.stem
            df = self.sketchy.predict_nanopore(
                fastq=fq,
                sketch=sketch,
                data=data,
                cores=cores,
                score=True,
                header=False,
                nreads=reads,
                top=1,
                out=None,
                tmp=self.outdir / f'{replicate_name}_tmp',
                sort_by='shared',
                quiet=True
            )

            # Bootstrap replicate IDs
            df['bootstrap'] = [replicate_name for _ in df.score]
            score_data.append(df)

            shutil.rmtree(self.outdir / f'{replicate_name}_tmp')

        return pandas.concat(score_data)

    def plot_bootstraps(
        self,
        bootstrap_data: Path,
        truth_data: Path = None,
        confidence=0.95,
        display=False
    ):

        df = pandas.read_csv(bootstrap_data, sep='\t', header=0, index_col=0)

        group = 0
        lineages, genotypes = [], []
        lineage_scores, genotype_scores = [], []
        failed_lineages, failed_genotypes = 0, 0
        for b, group in df.groupby('bootstrap'):
            lineage, genotype = self.get_thresholds(
                df=group,
                true_lineage=8,
                true_genotype="KL106-O2v2-299-0-0-0",
                true_susceptibility=None
            )

            lineage_score, genotype_score = self.get_line_data(
                df=group,
                true_lineage=8,
                true_genotype="KL106-O2v2-299-0-0-0"
            )

            lineage_scores.append(lineage_score)
            genotype_scores.append(genotype_score)

            lineages.append(lineage)
            genotypes.append(genotype)

            if lineage is None:
                failed_lineages += 1
            if genotype is None:
                failed_genotypes += 1

        if not lineages or not genotypes:
            raise RuntimeError(
                'Something went terribly wrong with the thresholds.'
            )

        bootstraps = dict(
            lineage=dict(
                data=lineages,
                ci=mean_confidence_interval(lineages, confidence=confidence),
                label="True Lineage"
            ),
            genotype=dict(
                data=genotypes,
                ci=mean_confidence_interval(genotypes, confidence=confidence),
                label="True Genotype"
            ),
            lineage_scores=dict(
                data=lineage_scores,
                label="Lineage Score Stability"
            ),
            genotype_scores=dict(
                data=genotype_scores,
                label="Genotype Score Stability"
            )
        )

        if failed_lineages > 0 or failed_genotypes > 0:
            print(
                f"Warning: there are {failed_lineages} failed lineages and "
                f"{failed_genotypes} failed genotypes. You should type a "
                f"larger number of reads (try > {2 * len(group)} reads)."
            )

        mu, lb, ub = bootstraps['lineage']['ci']
        print(
            f'{int(confidence * 100)}% CI lineage detection:'
            f' {int(lb)} - {int(mu)} - {int(ub)} reads'
        )

        mu, lb, ub = bootstraps['genotype']['ci']
        print(
            f'{int(confidence * 100)}% CI genotype detection:'
            f' {int(lb)} - {int(mu)} - {int(ub)} reads'
        )

        fig, (a1, a2) = plt.subplots(nrows=2, ncols=2)

        self.get_bootstrap_distplot(
            bootstraps, ax=a1[1]
        )
        self.get_bootstrap_ciplot(
            bootstraps, ax=a1[0]
        )
        self.get_confidence_interval_heatmap(
            bootstraps, ax=a2[0], confidence=confidence
        )

        self.get_bootstrap_lineplot(
            bootstraps, ax=a2[1]
        )

        if display:
            plt.show()
        else:
            plt.savefig('bootstrap_ci_line.pdf')
        plt.close()

    @staticmethod
    def get_line_data(df, true_lineage, true_genotype):

        lin = df.loc[df['primary_lineage'] == true_lineage]
        gen = df.loc[df['genotype'] == true_genotype]

        lineage, geno = [], []
        for bs, _ in enumerate(df.index):
            if bs in lin.index:
                lineage.append(lin.loc[bs, 'score'])
            else:
                lineage.append(0)

            if bs in gen.index:
                geno.append(gen.loc[bs, 'score'])
            else:
                geno.append(0)

        #  assert len(df) == len(line)

        return lineage, geno

    def get_bootstrap_lineplot(
            self,
            bootstraps: dict,
            ax=None
    ):

        self.make_line_plot(
            bootstraps, score_data='lineage_scores', palette='YlGnBu', ax=ax
        )
        genotype_plot = self.make_line_plot(
            bootstraps, score_data='genotype_scores', palette='YlOrBr', ax=ax
        )

        ax.set_xlabel('Reads')
        ax.set_ylabel('Mean preference score\nfor ST258')

        return genotype_plot

    @staticmethod
    def make_line_plot(
            bootstraps, score_data='lineage_scores', palette='YlGnBu',
            ax=None
    ):

        import operator
        from functools import reduce

        scores = reduce(
            operator.concat, bootstraps[score_data]['data']
        )

        bootstraps = len(
            bootstraps[score_data]['data']
        )

        reads = len(scores) // bootstraps
        read_vector = [i for i in range(reads)]

        read_score_vector = [read_vector for _ in range(bootstraps)]

        flattend_reads = reduce(
            operator.concat, read_score_vector
        )

        df = pandas.DataFrame(
            data={
                "scores": scores,
                "reads": flattend_reads
            }
        )

        sns.lineplot(
            data=df, x='reads', y='scores', lw=1.2,
            estimator='mean', n_boot=100, ci=None,
            color=sns.color_palette(palette, 6)[4],
            ax=ax
        )

        mean_bootstrap_scores = sns.lineplot(
            data=df, x='reads', y='scores', lw=1.2,
            estimator='std', n_boot=100, ci=None,
            color=sns.color_palette(palette, 6)[2],
            ax=ax
        )

        return mean_bootstrap_scores

    @staticmethod
    def get_bootstrap_ciplot(
            bootstraps: dict, barplot=True, hue=None, show_true=True,
            ax=None
    ):

        sns.set(style="white")

        df = pandas.DataFrame(
            data={
                'lineage': bootstraps['lineage']['ci'],
                'genotype':  bootstraps['genotype']['ci'],
            }
        )

        bs = pandas.DataFrame(
            data={
                'reads': bootstraps['lineage']['data'] +
                         bootstraps['genotype']['data'],
                'label': ['Lineage' for _ in bootstraps['lineage']['data']] +
                         ['Genotype' for _ in bootstraps['genotype']['data']],
            }
        )

        means = df.iloc[0, :]
        yerr = df.iloc[[1, 2], :]
        yerr = [
            np.absolute(
                (row - means)
            ) for i, row in yerr.iterrows()
        ]

        if barplot:
            if show_true:
                sns.stripplot(
                    data=pandas.DataFrame(data={
                        'reads': [10, 67],
                        'label': ['Lineage', 'Genotype']
                    }), hue=hue,
                    x='reads', y='label', jitter=False, size=5,
                    linewidth=1, palette="YlGnBu", alpha=.75, ax=ax
                )
            p = sns.barplot(
                data=bs, x='reads', y='label', xerr=yerr, hue=hue,
                linewidth=2.5, facecolor=(1, 1, 1, 0),
                errcolor=".2", edgecolor=".2", capsize=.05, ax=ax
            )

            ax.set_xlabel('Reads')
            ax.set_ylabel('')

            return p

        else:
            return plt.errorbar(
                x=df.columns.tolist(), y=means.tolist(), yerr=yerr, fmt='.k'
            )

    @staticmethod
    def get_confidence_interval_heatmap(
            bootstraps: dict,
            confidence: float = 0.95,
            ax=None
    ):
        sns.set(style="white")

        lci = list(bootstraps['lineage']['ci'])
        gci = list(bootstraps['genotype']['ci'])

        data = pandas.DataFrame(np.array(
            [[int(lci[1]), int(lci[0]), int(lci[2])],
            [int(gci[1]), int(gci[0]), int(gci[2])]]
        ), columns=['Lower', 'Mean', 'Upper'], index=['Lineage', 'Genotype'])

        p = sns.heatmap(
            data=data,
            cmap=sns.light_palette('#D3D3D3', n_colors=6),
            annot=True,
            linewidths=.5,
            fmt='.0f',
            cbar=False,
            ax=ax,
        )

        ax.set_xlabel(f"\n{int(confidence * 100)}% Confidence Interval (Reads)")
        ax.set_ylabel('')

        ax.tick_params(axis='y', labelrotation=0)

        return p


    @staticmethod
    def get_bootstrap_distplot(
        bootstraps: dict, boxplot=False, violinplot=True, hue=None, ax=None
    ):

        sns.set(style="white")

        df = pandas.DataFrame(
            data={
                'reads': bootstraps['lineage']['data'] +
                         bootstraps['genotype']['data'],
                'label': ['Lineage' for _ in bootstraps['lineage']['data']] +
                         ['Genotype' for _ in bootstraps['genotype']['data']],
            }
        )

        if boxplot:
            sns.boxplot(
                data=df,
                hue=hue,
                x='reads',
                y='label',
                palette="YlGnBu",
                ax=ax
            )
        if violinplot:
            sns.violinplot(
                data=df,
                hue=hue,
                x='reads',
                y='label',
                color=".8",
                ax=ax
            )

        p = sns.stripplot(
            data=df,
            hue=hue,
            x='reads',
            y='label',
            jitter=True,
            size=2.5,
            linewidth=1,
            palette="YlGnBu",
            alpha=.25,
            ax=ax
        )

        ax.set_xlabel('Reads')
        ax.set_ylabel('')

        return p

    @staticmethod
    def get_thresholds(
        df: pandas.DataFrame,
        true_lineage=None,
        true_genotype=None,
        true_susceptibility=None
    ):

        lin1 = df.loc[df['primary_lineage'] == true_lineage]
        try:
            first_lineage = lin1.index[0]
        except IndexError:
            first_lineage = None

        lin2 = df.loc[df['genotype'] == true_genotype]
        try:
            first_genotype = lin2.index[0]
        except IndexError:
            first_genotype = None

        return first_lineage, first_genotype







