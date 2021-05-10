__version__ = '0.5.0'

import os
import shutil
import logging
import pandas
import json
import seaborn as sns
import matplotlib.pyplot as plt

from numpy import nan
from pathlib import Path
from sketchy.utils import run_cmd, PoreLogger
from collections import OrderedDict
from colorama import Fore

RE = Fore.RESET
C = Fore.CYAN
G = Fore.GREEN
Y = Fore.YELLOW
R = Fore.RED


class SketchyDiagnostics(PoreLogger):

    def __init__(
        self,
        outdir: Path = Path("diagnostics"),
        mpl_backend: str = None,
        verbose: bool = True
    ):

        self.na = "#d3d3d3"
        self.outdir = outdir
        self.outdir.mkdir(parents=True, exist_ok=True)

        if mpl_backend:
            plt.switch_backend(mpl_backend)

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

    def plot_genotype_heatmap(self, nextflow: Path, subset_column: str, subset_values: str):

        """ Main access function for comparative feature heatmaps from Nextflow """

        nextflow_files = nextflow.glob("*.tsv")

        scale = 1.0

        for file in nextflow_files:
            nxf = pandas.read_csv(file, sep="\t", index_col=0, header=0)

            if subset_column:
                sv = [_.strip() for _ in subset_values.split(',')]
                nxf = nxf[nxf[subset_column].isin(sv)]

            mode = nxf['mode'].unique()[0]  # get the analysis mode from the results

            if mode == "stream":
                nxf = nxf.drop(columns="read")  # drop unused read column
            elif mode == "dist":
                nxf = nxf.drop(columns=["rank", "distance", "shared_hashes"])  # drop unused read column
            elif mode == "screen":
                nxf = nxf.drop(columns=["rank", "identity", "shared_hashes"])  # drop unused read column
            else:
                raise ValueError(f"Unsupported mode: {mode}")

            for db, db_data in nxf.groupby('db'):
                self.logger.info(f"Processing database: {db}")\

                nrows = len(db_data["read_limit"].unique())
                fig, axes = plt.subplots(
                    nrows=nrows, ncols=1, figsize=(
                        1 * 4 * 9 * scale, nrows * 2 * 9 * scale
                    )
                )

                for (i, (read_limit, predictions)) in enumerate(
                    db_data.groupby('read_limit')
                ):
                    _drop_for_labels = ['db', 'mode', 'read_limit']
                    if mode in ('dist', 'screen'):
                        _drop_for_labels.append('id')

                    _predictions = predictions.drop(columns=_drop_for_labels)

                    self.plot_comparative_heatmap(
                        values=None, annot=True, cbar=False,
                        labels=_predictions.sort_index(), palette="Set2_r",
                        title=f"\n{read_limit} Reads\n", ax=axes[i]
                    )

                plt.tight_layout()
                fig.savefig(f"{self.outdir / f'{mode}.{db}'}.png")

    def plot_sssh_diagostics(
        self, sssh_data: dict, plot_file: Path, plot_breakpoint: bool = False, color: str = "YlGnBu"
    ):

        # Something odd with colors, need reverse palettes:
        if not color.endswith('_r'):
            color += '_r'

        number_features = len(sssh_data)
        number_plots = 2

        fig, axes = plt.subplots(
            nrows=number_features, ncols=number_plots, figsize=(
                number_plots * 7, number_features * 4.5
            )
        )

        fig.subplots_adjust(hspace=0.8)

        for (i, (feature, data)) in enumerate(sssh_data.items()):

            feature_data = data["feature_data"]

            self.plot_sssh(
                feature_name=feature,
                feature_data=feature_data,
                top_feature_values=data['feature_values'],
                stability_breakpoint=data['stable_breakpoint'],
                color=color,
                plot_breakpoint=plot_breakpoint,
                ax=axes[i, 0]
            )

            self.plot_preference_score(
                feature_data=feature_data[feature_data['feature_rank'] == 0].reset_index(),
                ax=axes[i, 1]
            )

        plt.tight_layout()
        fig.savefig(f"{self.outdir / plot_file}")

    def plot_ssh_diagnostics(
        self, db: Path, ssh_file: Path, plot_file: Path, max_ranks: int = 5, color: str = "YlGnBu"
    ):

        ssh = self.read_ssh(file=ssh_file)
        genotypes = self.read_genotypes(db=db)

        ssh_features = ssh \
            .join(genotypes, how='inner') \
            .sort_values(['read', 'rank'])

        reads = len(ssh_features['read'].unique())
        ranks = len(ssh_features['rank'].unique())

        features = [_ for _ in genotypes.columns.tolist() if _ != 'idx']

        fig, ax = plt.subplots(
            nrows=len(features), ncols=1, figsize=(
                1 * 7, len(features) * 4.5
            )
        )

        fig.subplots_adjust(hspace=0.8)

        for (i, feature_name) in enumerate(features):
            self.plot_ssh(
                ssh_features=ssh_features,
                feature_name=feature_name,
                max_ranks=max_ranks,
                reads=reads,
                ranks=ranks,
                color=color,
                ax=ax[i]
            )

        plt.tight_layout()
        fig.savefig(f"{self.outdir / plot_file}")

    @staticmethod
    def read_genotypes(db: Path):

        genotype_file = db / f"{db.name}.tsv"

        if not genotype_file.exists():
            raise ValueError(f"Could not read genotype file from database: {genotype_file}")

        df = pandas.read_csv(genotype_file, sep="\t", header=0)
        df['idx'] = [i for i in df.iterrows()]  # always sorted by ascending index from database create

        if 'id' not in df.columns.tolist():
            raise ValueError("Genotype file must have a column with genome identifiers named 'id'")

        df = df.drop(columns='id')

        return df

    @staticmethod
    def plot_preference_score(
        feature_data: pandas.DataFrame,
        threshold: float = 0.6,
        ax: plt.axes = None
    ) -> None:

        p3 = sns.lineplot(
            data=feature_data, x='read', y='preference_score',
            ax=ax, color='#333333', ci=None, estimator=None
        )

        ax.axhline(
            y=threshold, linewidth=1, linestyle='--', color='black'
        )

        # Legend and labels
        p3.tick_params(labelsize=10)
        p3.set_xlabel('\nReads', fontsize=12)
        p3.set_ylabel('Preference score\n', fontsize=12)

    @staticmethod
    def plot_comparative_heatmap(
        values: pandas.DataFrame = None,
        palette: [str] or str = "YlGnBu",
        ax=None,
        fmt: str = ".3f",
        cbar: bool = True,
        annot: bool = True,
        labels: pandas.DataFrame = None,
        threshold: float = 0.,
        time: bool = False,
        title: str = "",
        evaluation: bool = False
    ):

        if not palette.endswith("_r"):
            palette += "_r"

        if values is None:
            if labels is None:
                raise ValueError("If no values supplied, a label matrix is required")
            values = labels.replace(labels, 1.)
        else:
            if all(values.isna().all().tolist()):  # if all values are NA
                values = values.fillna(0.)

        p1 = sns.heatmap(
            values,
            vmin=0 if evaluation else None,
            vmax=3 if evaluation else None,
            linewidths=5,
            cbar=cbar,
            ax=ax,
            annot=annot,
            fmt=fmt,
            cmap=palette,
            annot_kws={"size": 18 if time else 30, "weight": "bold"}
        )
        p1.tick_params(axis='both', which='major', labelsize=30, length=3, width=2)
        p1.tick_params(axis='x', rotation=90)
        p1.tick_params(axis='y', rotation=0)

        ax.set_title(title, fontdict={'fontsize': 38})
        ax.set_xlabel('')

        if threshold > 0.:
            for text in ax.texts:
                if float(text.get_text()) < threshold:
                    text.set_text("")

        if not time and labels is not None:
            label_vec = labels.stack().tolist()
            for i, text in enumerate(ax.texts):
                try:
                    # there is always categorical data never numeric floats
                    val = f"{float(label_vec[i]):.0f}"
                except ValueError:
                    val = label_vec[i].strip()

                    # Specific database rules for better visibility
                    if val.startswith('SCCmec'):
                        val = val.replace('SCCmec-', '')

                text.set_text(val)

        if time:
            label_vec = labels.stack().tolist()
            for i, text in enumerate(ax.texts):
                if isinstance(label_vec[i], str):
                    try:
                        text.set_text(
                            label_vec[i].split(" ")[1]
                        )  # date time
                    except KeyError:
                        text.set_text(
                            label_vec[i]
                        )  # time delta
                else:
                    text.set_text()

    def process_sssh(
        self, sssh_file: Path, stable: int = 100, mode: str = "last", max_ranks: int = 5
    ) -> dict:

        """ Process the raw output of the predict subcommand """

        sssh = self.read_sssh(sssh_file)

        sssh_data = {}
        for (i, (feature, feature_data)) in enumerate(sssh.groupby('feature')):  # each distinct feature is processed

            stable_breakpoint = self.compute_breakpoint(feature_data=feature_data, break_point=stable)

            feature_prediction, feature_values, preference_score = self.get_feature_prediction(
                feature_data, mode=mode, max_ranks=max_ranks
            )

            sssh_data[feature] = {
                'stable_breakpoint': stable_breakpoint,
                'feature_prediction': feature_prediction,
                'preference_score': preference_score,
                'feature_values': feature_values,
                'feature_data': feature_data,
                'max_ranks': max_ranks
            }

            feature_data.to_csv(self.outdir / f"{feature}.tsv", sep="\t", index=False, header=True)

        return sssh_data

    @staticmethod
    def plot_sssh(
        feature_name: str,
        feature_data: pandas.DataFrame,
        top_feature_values: list,
        stability_breakpoint: int = None,
        ax: plt.axes = None,
        color: str = 'YlGnBu',
        plot_breakpoint: bool = False,
    ) -> None:

        feature_data = feature_data.loc[
            feature_data['feature_value'].isin(top_feature_values), :
        ]

        feature_data = feature_data.assign(
            feature_value=feature_data['feature_value'].astype(str)
        )

        feature_values = feature_data.feature_value.unique()

        palette = sns.color_palette(
            color, n_colors=len(feature_values)
        )[:len(top_feature_values)]

        p2 = sns.lineplot(
            data=feature_data, x='read', y='sssh', hue='feature_value',
            ax=ax, ci=None, estimator=None, palette=palette,
            hue_order=[str(v) for v in top_feature_values]
        )

        if stability_breakpoint is not None and plot_breakpoint:
            ax.axvline(
                x=stability_breakpoint, linewidth=1, linestyle='-', color='black'
            )

        legend = ax.legend()
        legend.texts[0].set_text(feature_name)

        ax.set_title(feature_name, fontdict={'fontsize': 24})

        # Legend and labels
        p2.tick_params(labelsize=10)
        p2.set_xlabel('\nReads', fontsize=12)
        p2.set_ylabel('Sum of ranked sum of shared hashes\n', fontsize=12)

    @staticmethod
    def get_feature_prediction(
        feature_data: pandas.DataFrame, mode: str = "last", max_ranks: int = 5
    ):

        """ Gets the """

        if mode == "total":
            sorted_values = feature_data.groupby('feature_value') \
                .sum().sort_values(by='sssh', ascending=False)

            feature_prediction = sorted_values.iloc[0, :].name
            feature_values = sorted_values[:max_ranks].index.tolist()

            preference_score = float(
                feature_data['score'].median()
            )

        else:
            last_prediction = feature_data[feature_data["read"] == feature_data["read"].max()]

            feature_prediction = last_prediction.iloc[0, :].feature_value  # last top ranked feature value

            _feature_values = [feature_prediction] + \
                feature_data["feature_value"].value_counts()[:max_ranks].index.tolist()

            feature_values = []  # preserves hue order in diagnostic plots
            for v in _feature_values:
                if v not in feature_values:
                    feature_values.append(v)

            preference_score = last_prediction.iloc[0, :].preference_score

        return feature_prediction, feature_values, preference_score

    @staticmethod
    def compute_breakpoint(feature_data: pandas.DataFrame, break_point: int = None):

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

        if break_point:
            read_index = stable_point - break_point
            # Conditions: break point validations
            if last_stable_block_index == 0:  # none detected
                read_index = -1
        else:
            read_index = -1

        return read_index

    def plot_ssh(
        self,
        ssh_features: pandas.DataFrame,
        feature_name: str,
        reads: int,
        ranks: int,
        max_ranks: int = 5,
        ax: plt.axes = None,
        color: str = 'YlGnBu'
    ):

        if not color.endswith("_r"):
            color += "_r"

        top_values = ssh_features[feature_name].value_counts()[:max_ranks].index.tolist()

        palette = sns.color_palette(
            color, max_ranks
        )[:len(top_values)]

        df = ssh_features.copy().reset_index()
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

        p1.set_facecolor(self.na)

        xticks, yticks = self._get_ticks(reads, ranks)

        p1.set_xticks(xticks)
        p1.set_xticklabels(xticks, rotation='vertical')

        p1.set_xlabel('\nReads', fontsize=10)
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

    @staticmethod
    def _get_ticks(reads, ranks):

        if reads <= 200:
            xticks = [i for i in range(0, reads + 1, 20)]
        elif 200 < reads <= 500:
            xticks = [i for i in range(0, reads + 1, 50)]
        elif 500 < reads <= 1500:
            xticks = [i for i in range(0, reads + 1, 200)]
        elif 1500 < reads <= 5000:
            xticks = [i for i in range(0, reads + 1, 1000)]
        elif 5000 < reads <= 15000:
            xticks = [i for i in range(0, reads + 1, 3000)]
        else:
            xticks = [i for i in range(0, reads + 1, 5000)]

        if ranks > 10:
            yticks = [i for i in range(0, ranks + 1, 10)]
        else:
            yticks = [i for i in range(0, ranks + 1)]

        return xticks, yticks

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
                'stability',         # stability breakpoint ...
                'preference_score',  # ... and preference score are for the top prediction
            ],
            dtype={
                'read': int,
                'feature': str,
                'feature_value': str,
                'feature_rank': int,
                'sssh': int,
                'stability': int,
                'preference_score': float
            }
        )

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

    def match_reference(self, nextflow, reference):
        """ Match predictions from collected Nextflow results to reference table """

        ref = pandas.read_csv(reference, sep="\t", header=0, index_col=0)

        unique_references = ref.index.unique().tolist()

        # Excluding reference columns
        to_exclude = []
        for col in ref.columns.tolist():
            column_values = ref[col].tolist()
            if all([v == "-" for v in column_values]):
                to_exclude.append(col)
                print(f'Excluding column: {col}')

        ref = ref.drop(columns=to_exclude)

        methods_summary = []
        for collected in nextflow.glob("*.tsv"):
            method = collected.stem
            data = pandas.read_csv(collected, sep="\t", header=0, index_col=0)

            if len(data['replicate'].unique()) > 1:
                is_bootstrapped = True
            else:
                is_bootstrapped = False

            print(f"Bootstrap replicates: {is_bootstrapped}")

            # Excluding isolates:
            unique_samples = data.index.unique().tolist()
            not_in_ref = list(set(unique_samples).difference(unique_references))
            print(
                f"Excluding samples ({method}) not in reference and data:"
                f" {', '.join(not_in_ref if not_in_ref else '-')}"
            )

            data = data[~data.index.isin(not_in_ref)]
            data = data[ref.columns.tolist() + ["db", "read_limit", "replicate"]]

            summary = []
            comparisons = {}
            for db, db_data in data.groupby("db"):
                for read_limit, read_data in db_data.groupby("read_limit"):
                    for _, row in read_data.iterrows():
                        sample = row.name
                        replicate = row['replicate']

                        row = row.drop(labels=['db', 'read_limit', 'replicate'])

                        sample_ref = ref.loc[sample, :]
                        comparison = pandas.DataFrame(
                            [row, sample_ref],
                            index=[f"call", f"reference"]
                        ).T

                        comparison["match"] = comparison["call"] == comparison["reference"]

                        true_calls = sum([1 for v in comparison['match'] if v == True])
                        total_calls = len(comparison)
                        true_percent = round((true_calls/total_calls)*100, 2)

                        genotypes = comparison.index.tolist()

                        if 'st' in genotypes:
                            true_st = comparison.loc['st', 'match']
                        elif 'mlst' in genotypes:
                            true_st = comparison.loc['mlst', 'match']
                        elif 'lineage' in genotypes:
                            true_st = comparison.loc['lineage', 'match']
                        else:
                            true_st = nan

                        summary.append(
                            [sample, db, read_limit, true_calls, total_calls, true_percent, true_st, replicate]
                        )
                        comparisons[f"{db}_{read_limit}_{sample}_{replicate}"] = comparison

            summary_df = pandas.DataFrame(
                summary, columns=['sample', 'db', 'read_limit', 'true_calls', 'total_calls', 'true_percent', 'true_st', 'replicate']
            ).sort_values(['db', 'sample', 'read_limit'])

            summary_df['method'] = [method for _ in summary_df.iterrows()]

            methods_summary.append(summary_df)

        df = pandas.concat(methods_summary).reset_index(drop=True)

        for db, db_data in df.groupby("db"):
            fig, axes = plt.subplots(
                nrows=1, ncols=1, figsize=(14, 10)
            )

            sns.barplot(data=db_data, x="read_limit", y="true_percent", hue="method", ax=axes, palette="colorblind")

            plt.tight_layout()
            fig.savefig(f"{db}.summary.png")

        print(df)

class SketchyDatabase(PoreLogger):

    def __init__(
        self,
        db_path: Path = None,
        sketch_file: Path = None,
        genotype_file: Path = None,
        lineage_column: str = 'mlst',
        verbose: bool = True
    ):

        PoreLogger.__init__(
            self, level=logging.INFO if verbose else logging.ERROR
        )

        self.db_path = db_path
        self.sketch_file = sketch_file
        self.genotype_file = genotype_file

        if db_path:
            db_name = db_path.name
            if not db_path.exists():
                # Try with env variable SKETCHY_PATH
                try:
                    self.db_path = os.environ['SKETCHY_PATH'] / db_name
                except KeyError:
                    self.logger.error("Could not find database path")
                    exit(1)

            self.sketch_file, self.genotype_file, self.genotype_index, self.genotype_key = self.get_sketch_files()

        if self.genotype_file:
            self.genotypes = pandas.read_csv(
                self.genotype_file, sep='\t', header=0
            )
        else:
            self.genotypes = None

        self.lineage_column = lineage_column

    def create_database(
        self, outdir: Path = "sketchy-db", drop: str = None, numeric: bool = False
    ):

        self.logger.info(f"Database sketch: {self.sketch_file}")
        self.logger.info(f"Genotype file: {self.genotype_file}")
        self.logger.info(f"Database directory: {outdir}")
        self.logger.info(f"Drop columns: {drop}")
        self.logger.info(f"Include numeric columns: {drop}")

        outdir.mkdir(parents=True, exist_ok=True)

        id_column: str = 'id'

        # Drop columns from genotypes

        if drop:
            genotypes = self.genotypes.drop(columns=drop.split(','))
        else:
            genotypes = self.genotypes.copy()

        if 'id' not in genotypes.columns.tolist():
            raise ValueError("Genotype data file must have a column with genome identifiers named 'id'")

        # Read sketch index and extract genome identifiers:

        sketch_info = self.get_sketch_info()

        _isolates_in_sketch = len(sketch_info)
        _isolates_in_genotypes = len(genotypes)

        if _isolates_in_sketch != _isolates_in_genotypes:
            self.logger.error(
                f"Number of isolates in sketch ({_isolates_in_sketch}) "
                f"does not match genotypes ({_isolates_in_genotypes})"
            )
            exit(1)

        _names_in_sketch = sketch_info['id'].tolist()
        _names_in_genotypes = genotypes[id_column].tolist()

        if len(set(_names_in_sketch)) != len(_names_in_sketch):
            self.logger.error("Duplicate identifiers in sketch file! Please replace before proceeding")
            exit(1)

        if len(set(_names_in_genotypes)) != len(_names_in_genotypes):
            self.logger.error("Duplicate identifiers in genotype file! Please replace before proceeding")
            exit(1)

        if set(_names_in_sketch) != set(_names_in_genotypes):
            self.logger.error("Genotype identifiers do not match identifiers in sketch (stem of assembly names)")
            exit(1)

        indexed_genotypes = genotypes.merge(
            sketch_info, left_on=id_column, right_on="id", how='inner'
        )

        genotypes_reference = indexed_genotypes.sort_values('idx').set_index('idx')

        self.logger.info(f"Transforming genotype columns for numeric index for Sketchy")

        genotype_index, genotype_keys = self.transform_columns(
            genotypes=indexed_genotypes.copy(), numeric=numeric
        )

        genotype_index['idx'] = indexed_genotypes['idx']
        genotype_index = genotype_index.sort_values('idx').set_index('idx')

        _path = outdir / outdir.name

        self.logger.info(f"Writing database files to: {outdir}")

        genotype_index.to_csv(_path.with_suffix('.idx'), index=False, sep="\t", header=False)
        genotypes_reference.to_csv(_path.with_suffix('.tsv'), index=False, sep="\t", header=True)

        with _path.with_suffix('.key').open('w') as fout:
            json.dump(genotype_keys, fout, sort_keys=False)

        shutil.copyfile(self.sketch_file, str(_path.with_suffix('.msh')))

    def transform_columns(self, genotypes: pandas.DataFrame, numeric: bool = True):

        """ Transform into categorical numeric data for evaluation in Rust """

        dtypes = ['category', 'bool', 'object']
        if numeric:  # also transform numeric values of type int64 (to categorical)
            dtypes += ['int64']

        transform = genotypes.select_dtypes(dtypes).columns

        transform = [_ for _ in transform if _ not in ('id', 'idx')]

        genotypes.drop(columns=[
            c for c in genotypes.columns if c not in transform
        ], inplace=True)

        feature_keys = OrderedDict()
        for (i, (name, column_data)) in enumerate(
            genotypes.iteritems()
        ):
            self.logger.info(f"Processing genotype: {name}")
            feature_keys[i] = {
                'name': name,
                'values': [str(v) for v in column_data.astype('category').cat.categories.tolist()]
            }

        genotypes[transform] = genotypes[transform].apply(
            lambda x: x.astype('category').cat.codes
        )

        return genotypes, feature_keys

    def get_sketch_files(self):

        db_name = self.db_path.name

        sketch_file = self.db_path / (db_name + ".msh")
        genotype_file = self.db_path / (db_name + ".tsv")
        genotype_index = self.db_path / (db_name + ".idx")
        genotype_key = self.db_path / (db_name + ".key")

        for f in (sketch_file, genotype_file, genotype_index, genotype_key):
            if not f.exists():
                self.logger.error(f"Database file does not exist: {f}")

        return sketch_file, genotype_file, genotype_index, genotype_key

    def get_sketch_info(self) -> pandas.DataFrame:

        run_cmd(f'mash info -t {self.sketch_file} > info.tmp', shell=True)

        converters = {'fname': lambda x: Path(x).stem}
        mash_info = pandas.read_csv(
            f'info.tmp',
            sep='\t',
            header=None,
            skiprows=1,
            index_col=0,
            engine='c',
            usecols=[2],
            names=['fname'],
            converters=converters,
        )
        os.remove('info.tmp')

        mash_info['idx'] = [i for i in range(0, len(mash_info))]
        mash_info['id'] = mash_info.index.tolist()

        return mash_info

    def write(self, file: Path, idx: bool = True, header: bool = True):

        self.genotypes.sort_values('idx', ascending=True).to_csv(
            file, sep='\t', header=header, index=idx
        )

    def has_lineage(self, lineage: str) -> bool:

        return lineage in self.genotypes[self.lineage_column].values

    def get_lineage(self, lineage: str) -> pandas.DataFrame:

        """" Get a subset of the genotype index for the given db_lineage """

        if self.has_lineage(lineage):
            return self.genotypes[self.genotypes[self.lineage_column] == lineage]
        else:
            raise ValueError(f'Could not detect lineage in index: {lineage}')

    def get_key_index(
        self, lineage: str, key_file: Path = None,
    ) -> pandas.DataFrame:

        """ Access feature data by legacy key index file from Pathfinder Survey """

        df = self.get_lineage(lineage=lineage)

        key = pandas.read_csv(
            key_file, sep="\t", usecols=['uuid', 'id'],
            header=0, index_col=None
        )

        keyed = df.merge(key, how='inner', on='uuid')
        keyed.index = df.index

        return keyed

    def drop_columns(self, columns: list = None):

        for column in columns:
            try:
                self.genotypes.drop(columns=column, inplace=True)
            except KeyError:
                print(f'Could not find column: {column} - skipping')

        return self.genotypes

    def get_summary(self, lineage: str) -> pandas.DataFrame or None:

        """ Generate a genotype summary for the given db_lineage """

        ignore_columns = ['id']

        if self.has_lineage(lineage):

            # Get frequency summary of create values

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
            self.logger.error(f'Could not detect lineage in index: {lineage}')
            exit(1)

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
