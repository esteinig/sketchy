__version__ = '0.4.4'

import logging
import pandas
import json

import seaborn as sns

import matplotlib.pyplot as plt

from numpy import nan
from numpy import reshape
from pathlib import Path
from sketchy.utils import run_cmd, is_fastq, PoreLogger
from collections import OrderedDict
from colorama import Fore

RE = Fore.RESET
C = Fore.CYAN
G = Fore.GREEN
Y = Fore.YELLOW
R = Fore.RED


class SketchyWrapper(PoreLogger):

    """ Python wrapper for the Rust pipeline """

    def __init__(
        self,
        fastx: Path,
        sketch: Path,
        prefix: str = 'sketchy',
        outdir: Path = Path('sketchy_out'),
        verbose: bool = False,
    ):

        PoreLogger.__init__(
            self,
            level=logging.INFO if verbose else logging.ERROR,
            name='Compute'
        )

        self.fastx = fastx
        self.sketch = sketch
        self.prefix = prefix
        self.outdir = outdir
        self.verbose = verbose

        self.logger.info(f'Sketchy wrapper v{__version__}')
        self.logger.info(f'Prefix: {prefix}')
        self.logger.info(f'Fastq file: {fastx.absolute()}')
        self.logger.info(f'Output directory: {outdir.absolute()}')

        self.outdir.mkdir(exist_ok=True, parents=True)

    def run(
        self,
        ranks: int = 10,
        limit: int = None,
        stable: int = 1000,
        threads: int = 4,
        plot: bool = True,
        palette: str = 'YlGnBu',
        image_format: str = 'pdf',
        mpl_backend: str = ""
    ) -> None:

        sketch, features, keys = self.get_sketch_files()

        if limit is not None:
            limit_pipe = f'| head -{limit* 4 if is_fastq(self.fastx) else 2}'
        else:
            limit_pipe = ''

        self.logger.info(f'Sketch database: {sketch}')
        self.logger.info(f'Consensus ranks: {ranks}')
        self.logger.info(f'Read limit: {"all" if limit is None else limit}')
        self.logger.info(f'Stability breakpoint: {stable}')
        self.logger.info(f'Threads for Mash: {threads}')

        command_ssh = f'cat {self.fastx} {limit_pipe}' \
            f' | sketchy-rs compute -r {ranks} -s {sketch}' \
            f' -t {threads} -p {1 if self.verbose else 0}' \
            f' > {self.outdir / self.prefix}.ssh.tsv'

        command_sssh = f'cat {self.outdir / self.prefix}.ssh.tsv' \
            f' | sketchy-rs evaluate -f {features} -s {stable}' \
            f' > {self.outdir / self.prefix}.sssh.tsv'

        self.logger.info('Computing sum of shared hashes...')
        run_cmd(command_ssh, shell=True)

        self.logger.info('Evaluating sum of shared hashes...')
        run_cmd(command_sssh, shell=True)

        eve = Evaluation(
            sssh=Path(f"{self.outdir / self.prefix}.sssh.tsv"),
            ssh=Path(f'{self.outdir / self.prefix}.ssh.tsv'),
            index=features,
            key=keys,
            stable=stable,
            verbose=self.verbose
        )

        eve.plot_feature_evaluations(
            plot_file=Path(f'{self.outdir / self.prefix}.{image_format}'),
            break_file=Path(f'{self.outdir / self.prefix}.data.tsv'),
            color=palette, break_point=True, plot=plot, mpl_backend=mpl_backend
        )

    def get_sketch_files(self):

        """ Find prefixed sketch collection or use local cache """

        mash = Path(str(self.sketch) + '.msh')
        feature = Path(str(self.sketch) + '.tsv')
        key = Path(str(self.sketch) + '.json')

        for f in (mash, feature, key):
            if not f.exists():
                raise ValueError(f'Could not detect sketch file: {f}')

        return mash, feature, key

    def check_rust_dependencies(self):

        """ Check dependency versions for Rust pipeline """

        try:
            output = run_cmd('rustc --version')
            rustc_version = output.decode('utf-8').split()[1].strip()
            self.logger.info(f'Rustc version: {rustc_version}')
        except FileNotFoundError:
            self.logger.info('Failed to run Sketchy: no `rustc` in $PATH.')
            exit(1)
        except KeyError:
            self.logger.info('Failed to parse version of: rustc')
            exit(1)

        try:
            output = run_cmd('mash --version')
            mash_version = output.decode('utf-8').strip()
            self.logger.info(f'Mash version: {mash_version}')
        except FileNotFoundError:
            self.logger.info('Failed to run Sketchy: no `mash` in $PATH.')
            exit(1)
        except KeyError:
            self.logger.info('Failed to parse version of: mash')
            exit(1)

        try:
            output = run_cmd('sketchy-rs --version')
            sketchyrs_version = output.decode('utf-8').split()[1].strip()
            self.logger.info(f'Sketchy Rust version: {sketchyrs_version}')
        except FileNotFoundError:
            self.logger.info('Failed to run Sketchy: no `sketchy-rs` in $PATH.')
            exit(1)
        except KeyError:
            self.logger.info('Failed to parse version of: sketchy-rs')
            exit(1)


class Evaluation(PoreLogger):

    """ Evaluations and plotting from Rust libraries """

    def __init__(
        self,
        sssh: Path,
        index: Path,
        key: Path,
        stable: int = None,
        ssh: Path = None,
        verbose: bool = False
    ):

        PoreLogger.__init__(self, name="Evaluate")

        if verbose:
            self.logger.setLevel(level=logging.INFO)

        self.top_feature_values = 5
        self.read_limit = 1000
        self.preference_threshold = 0.6
        self.na_color = 'darkgray'

        self.stable = stable

        self.logger.info(f"Loading data for evaluations from Sketchy Rust")
        self.logger.info(f"Ranked sum of shared hashes: {sssh}")
        self.logger.info(f"Sum of shared hashes: {ssh}")
        self.logger.info(f"Genotype feature index: {index}")
        self.logger.info(f"Genotype feature key: {key}")

        self.feature_key = self.read_feature_key(file=key)  # key to headers and categories
        self.feature_index, self.feature_data = self.read_feature_index(file=index)

        self.ssh = self.read_ssh(file=ssh)
        self.sssh = self.read_sssh(file=sssh)

        self.features = self.feature_index.columns.tolist()

        if self.ssh is not None:
            # Merge ssh and feature index for heatmap
            self.ssh_features = self.ssh \
                .join(self.feature_data, how='inner') \
                .sort_values(['read', 'rank'])

            self.reads = len(
                self.ssh_features['read'].unique()
            )

            self.ranks = len(
                self.ssh_features['rank'].unique()
            )

    def plot_feature_evaluations(
        self,
        plot_file: Path,
        break_file: Path,
        plot: bool = True,
        color: str = "YlGnBu",
        mpl_backend: str = "",
        break_point: bool = False
    ):

        self.logger.info(f"Compute and plot feature evaluations")
        self.logger.info(f"Plot break point: {break_point}")
        self.logger.info(f"Color palette: {color}")
        self.logger.info(f"Output predictions to file: {break_file}")

        # Something odd with colors, need reverse palettes:
        if not color.endswith('_r'):
            color += '_r'

        if plot:
            self.logger.info(f"Output plot to file: {plot_file}")
            number_features = len(self.features)
            number_plots = 3 if self.ssh is not None else 2
            fig, axes = plt.subplots(
                nrows=number_features, ncols=number_plots, figsize=(
                    number_plots * 7, number_features * 4.5
                )
            )

            if axes.ndim == 1:
                axes = reshape(
                    axes, (-1, 2)
                )

            fig.subplots_adjust(hspace=0.8)

        data = {}
        for (i, (feature, feature_data)) in enumerate(
            self.sssh.groupby('feature')
        ):

            feature_data, feature_name = self.translate_feature_data(
                df=feature_data, feature=str(feature)
            )

            top_prediction, top_values = self.get_top_feature_data(
                feature_data=feature_data
            )

            stability_breakpoint = self.compute_breakpoint(feature_data)

            median_preference = float(
                feature_data['score'].median()
            )

            data[feature_name] = {
                'stability': stability_breakpoint,
                'prediction': top_prediction,
                'preference': median_preference
            }

            if plot:

                if mpl_backend:
                    plt.switch_backend(mpl_backend)

                if self.ssh is not None:
                    self.plot_heatmap(
                        feature_name=feature_name,
                        top_values=top_values,
                        color=color,
                        ax=axes[i, 0]
                    )

                self.plot_sssh(
                    feature_name=feature_name,
                    feature_data=feature_data,
                    top_feature_values=top_values,
                    stability_breakpoint=stability_breakpoint,
                    color=color,
                    break_point=break_point,
                    ax=axes[i, 0 if self.ssh is None else 1]
                )

                single_score_data = feature_data[
                    feature_data['feature_rank'] == 0
                ].reset_index()

                self.plot_preference_score(
                    feature_data=single_score_data,
                    ax=axes[i, 1 if self.ssh is None else 2]
                )

                self.logger.info(f"Constructed plots for feature: {feature_name}")

        break_data = pandas.DataFrame(data).T
        break_data = break_data[['prediction', 'stability', 'preference']]

        break_data.stability = break_data.stability.astype(int)+1
        break_data.to_csv(break_file, sep='\t', index=True, index_label="feature")

        if plot:
            plt.tight_layout()
            fig.savefig(plot_file)
            self.logger.info(f'Saved evaluation plot to: {plot_file}')

        self.logger.info(f'Saved predictions to: {break_file}')

    def get_break_time_data(self, file: Path, data: pandas.DataFrame):

        self.logger.info('Determining breakpoint time from reads.')
        self.logger.info(f'Read file: {file}')

    def get_top_feature_data(self, feature_data: pandas.DataFrame):

        sorted_values = feature_data.groupby('feature_value') \
            .sum().sort_values(by='sssh', ascending=False)

        top_feature_prediction = sorted_values.iloc[0, :].name
        top_feature_values = sorted_values[:self.top_feature_values].index.tolist()

        return top_feature_prediction, top_feature_values

    @staticmethod
    def read_feature_key(file: Path):

        with file.open('r') as key_file:
            return json.load(key_file)

    def read_feature_index(self, file: Path):

        df = pandas.read_csv(
            file,
            sep='\t',
            header=None,
            index_col=False,  # ordered by sequential idx
            dtype='Int64'
        )

        feature_index = df.copy()

        column_names = []
        for name, column in df.iteritems():

            try:
                column_name = self.feature_key[str(name)]['name']
            except KeyError:
                raise KeyError(f'Could not get column name from feature {name}')

            column = column.astype('category')

            try:
                if -1 in column.cat.categories:
                    column.cat.categories = ['-1'] + self.feature_key[str(name)]['values']
                else:
                    column.cat.categories = self.feature_key[str(name)]['values']
            except KeyError:
                raise KeyError(
                    f'Could not find {name} or associated values in feature index key.'
                )

            df[name] = column
            column_names.append(column_name)

        df.index.name, feature_data = 'idx', df.copy()

        feature_data.columns = column_names
        feature_index.columns = column_names

        return feature_index, feature_data

    def compute_breakpoint(self, feature_data: pandas.DataFrame):

        top_predictions = feature_data[feature_data['feature_rank'] == 0]
        reverse_stability = top_predictions.stability.values[::-1]

        last_stable_block_index = 0
        for stable in reverse_stability:
            if stable == 0:
                break
            else:
                last_stable_block_index += 1

        # total reads - block length
        stable_point = len(reverse_stability) - last_stable_block_index+1

        if self.stable:
            read_index = stable_point - self.stable
            # Conditions: break point validations
            if last_stable_block_index == 0:  # none detected
                read_index = -1
        else:
            read_index = -1

        return read_index

    @staticmethod
    def read_ssh(file: Path = None):

        if file is None:
            return None
        else:
            df = pandas.read_csv(
                file, sep='\t', index_col=False, header=None,
                names=['idx', 'ssh', 'rank', 'read'], dtype={
                    'idx': int, 'ssh': int, 'rank': int, 'read': int
                }
            )

            return df.set_index('idx')

    @staticmethod
    def read_sssh(file: Path):

        return pandas.read_csv(
            file,
            sep='\t',
            header=None,
            index_col=False,  # reads as first column rather than index
            names=[
                'read',
                'feature',
                'feature_value',
                'feature_rank',
                'sssh',
                'stability',
                'score',
            ],
            dtype={
                'read': int,
                'feature': int,
                'feature_value': int,
                'feature_rank': int,
                'sssh': int,
                'stability': int,
                'score': float
            }
        )

    def plot_preference_score(
        self,
        feature_data: pandas.DataFrame,
        ax: plt.axes = None
    ) -> None:

        p3 = sns.lineplot(
            data=feature_data, x='read', y='score',
            ax=ax, color='#333333', ci=None, estimator=None
        )

        ax.axhline(
            y=self.preference_threshold, linewidth=1, linestyle='--', color='black'
        )

        # Legend and labels
        p3.tick_params(labelsize=6)
        p3.set_xlabel('\nReads', fontsize=9)
        p3.set_ylabel('Preference score\n', fontsize=9)

    def plot_sssh(
        self,
        feature_name: str,
        feature_data: pandas.DataFrame,
        top_feature_values: list,
        stability_breakpoint: int = None,
        ax: plt.axes = None,
        color: str = 'YlGnBu',
        break_point: bool = False,
    ) -> None:

        feature_data = feature_data.loc[
            feature_data['feature_value'].isin(top_feature_values), :
        ]

        feature_data = feature_data.assign(
            feature_value=feature_data['feature_value'].astype(str)
        )

        feature_values = feature_data.feature_value.unique()
        palette = sns.color_palette(
            color, n_colors=self.top_feature_values
        )[:len(feature_values)]

        p2 = sns.lineplot(
            data=feature_data, x='read', y='sssh', hue='feature_value',
            ax=ax, ci=None, estimator=None, palette=palette,
            hue_order=[str(v) for v in top_feature_values]
        )

        if stability_breakpoint is not None and break_point:
            ax.axvline(
                x=stability_breakpoint, linewidth=1, linestyle='-', color='black'
            )

        legend = ax.legend()
        legend.texts[0].set_text(feature_name)

        # Legend and labels
        p2.tick_params(labelsize=6)
        p2.set_xlabel('\nReads', fontsize=9)
        p2.set_ylabel('Sum of ranked sum of shared hashes\n', fontsize=9)

    def plot_heatmap(
        self,
        feature_name: str,
        top_values: list,
        ax: plt.axes = None,
        color: str = 'YlGnBu'
    ):

        palette = sns.color_palette(
            color, self.top_feature_values
        )[:len(top_values)]

        df = self.ssh_features.copy().reset_index()
        df = df[['rank', 'read', feature_name]]

        # Replace numeric for top features
        df = self._replace_heatmap_column(
            df=df, feature_name=feature_name, top_values=top_values
        )

        hm = df.pivot('rank', 'read', feature_name)

        hm = hm[hm.columns].astype(float)

        p1 = sns.heatmap(
            hm.iloc[:, :], linewidths=0, cbar=False, ax=ax, cmap=palette
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

    @staticmethod
    def _replace_heatmap_column(
        df: pandas.DataFrame,
        feature_name: str,
        top_values: list
    ):

        values = df[feature_name].unique()

        to_replace = {
            f: (top_values.index(f) if f in top_values else nan) for f in values
        }

        df[feature_name].replace(to_replace=to_replace, inplace=True)

        return df

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

    def translate_feature(self, feature_name, feature_data):

        pass

    def translate_feature_data(
        self, df: pandas.DataFrame, feature: str
    ):

        """ Category replacement should be fast """

        try:
            feature_dict = self.feature_key[feature]
        except KeyError:
            raise KeyError(
                f'Could not find translation data for feature {feature}'
            )

        try:
            feature_name = feature_dict['name']
            feature_values = {
                i: val for i, val in enumerate(
                    feature_dict['values']
                )
            }
        except KeyError:
            raise KeyError(
                f"Could not find name or value keys for feature: {feature}"
            )

        # Feature column
        df = df.assign(
            feature=df.feature.astype('category')
        )
        df.feature.cat.categories = [feature_name]

        # Feature value column
        df = df.assign(
            feature_value=df.feature_value.astype('category')
        )
        df.feature_value.cat.categories = [
            feature_values[cat] for cat in df.feature_value.cat.categories
        ]

        return df, feature_name


class LineageIndex(PoreLogger):

    def __init__(
        self,
        index_file: Path = None,
        lineage_column: str = 'mlst',
        verbose: bool = True
    ):

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

        if index_file:
            self.index = pandas.read_csv(
                index_file, sep='\t', header=0
            )
            if 'idx' not in self.index.columns:
                self.logger.info('Adding Mash index column "idx" to genotype index')
                self.index.index = [i for i in range(
                    len(self.index)
                )]
                self.index.index.name = 'idx'
            else:
                self.logger.info('Using column "idx" as Mash index column in genotype index')
                self.index.index = self.index.set_index('idx')
        else:
            self.index = None

        self.lineage_column = lineage_column

    def write(self, file: Path, idx: bool = True, header: bool = True):

        self.index.sort_values('idx', ascending=True).to_csv(
            file, sep='\t', header=header, index=idx
        )

    def has_lineage(self, lineage: str) -> bool:

        return lineage in self.index[self.lineage_column].values

    def get_lineage(self, lineage: str) -> pandas.DataFrame:

        """" Get a subset of the genotype index for the given lineage """

        if self.has_lineage(lineage):
            return self.index[self.index[self.lineage_column] == lineage]
        else:
            raise ValueError(f'Could not detect lineage in index: {lineage}')

    def get_key_index(
        self, lineage: str, key_file: Path = None,
    ) -> pandas.DataFrame:

        """ Access lineage data by legacy key index file from Pathfinder Survey """

        df = self.get_lineage(lineage=lineage)

        key = pandas.read_csv(
            key_file, sep="\t", usecols=['uuid', 'id'],
            header=0, index_col=None
        )

        keyed = df.merge(key, how='inner', on='uuid')
        keyed.index = df.index

        return keyed

    def drop_columns(self, columns: list = None):

        if columns is None:
            columns = ['uuid', 'id']

        for column in columns:
            try:
                self.index.drop(columns=column, inplace=True)
            except KeyError:
                print(f'Could not find column: {column} - skipping')

        return self.index

    def prepare_columns(self, integers: bool = True):

        """ Transform into categorical numeric data for evaluation in Rust """

        dtypes = ['category', 'bool', 'object']
        if integers:
            dtypes += ['int64']

        transform = self.index.select_dtypes(dtypes).columns

        self.index.drop(columns=[
            c for c in self.index.columns if c not in transform
        ], inplace=True)

        feature_keys = OrderedDict()
        for (i, (name, column_data)) in enumerate(
            self.index.iteritems()
        ):
            feature_keys[i] = {
                'name': name,
                'values': column_data.astype('category').cat.categories.tolist()
            }

        self.index[transform] = self.index[transform].apply(
            lambda x: x.astype('category').cat.codes
        )

        return self.index, feature_keys

    def get_summary(self, lineage: str) -> pandas.DataFrame or None:

        """ Generate a genotype summary for the given lineage """

        ignore_columns = ['uuid', 'idx', 'id', 'plasmid', 'resistance']  # TODO check these

        if self.has_lineage(lineage):

            # Get frequency summary of feature values

            df = self.get_lineage(lineage)

            feature_data = []
            for column in df.columns:
                if column not in ignore_columns:
                    feature_counts = df[column].value_counts()
                    feature_str, feature_df = self._build_summary_str(feature_counts)
                    feature_data.append(feature_df)
                    print(feature_str)

            return pandas.concat(feature_data)

        else:
            raise ValueError(f'Could not detect lineage in index: {lineage}')

    @staticmethod
    def _build_summary_str(feature_counts: pandas.Series):

        """ Build a summary string from value counts of a feature """

        feature_str = f"\n{G}{feature_counts.name.upper()}{RE}\n\n"

        total = sum(feature_counts)
        summary_data = []
        for idx, count in feature_counts.items():
            idx = idx.strip()

            percent = round(
                (count/total)*100, ndigits=2
            )

            if len(idx) >= 16:
                _idx = idx[:12] + "..."
            else:
                _idx = idx

            feature_str += f"{_idx:<16}{count:<8}{percent:8}%\n"
            summary_data.append(
                [feature_counts.name, idx, count, percent]
            )

        df = pandas.DataFrame(
            summary_data, columns=['feature', 'genotype', 'count', 'percent']
        )

        return feature_str, df
