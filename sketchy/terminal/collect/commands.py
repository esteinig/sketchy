import click
import pandas

from pathlib import Path
from matplotlib import pyplot as plt
from numpy import reshape
import seaborn as sns


@click.command()
@click.option(
    '--directory', '-d', type=Path, help='Path to directory to collect {prefix}.data.tsv'
)
@click.option(
    '--nextflow', '-n', is_flag=True, help='Predictions are from Sketchy Nextflow'
)
@click.option(
    '--prefix', '-p', type=str, default="summary", help='Prefix for summary files [summary]'
)
@click.option(
    '--subset', '-s', type=str, default=None,
    help='When using Nextflow use a subset string for specific configurations'
         'of ranks & reads: 10,1000 - or a sample prefix: isolate1  [None]'
)
@click.option(
    '--reference', '-r', type=str, default=None,
    help='Genotype matrix in same format as output containing feature truths'
)
@click.option(
    '--heatmap', '-m', is_flag=True,
    help="Visualize results as heatmap"
)
@click.option(
    '--time', '-st', is_flag=True,
    help="Parse the time enhanced output files for Nextflow"
)
@click.option(
    '--threshold', '-t', type=float, default=0.6,
    help="Apply threshold value to median preference score summary; values below are set to 0 [0.6] "
)
def collect(directory, nextflow, prefix, subset, heatmap, threshold, time, reference):

    """ Collect predictions and summarize results """

    if subset:
        subset = subset.split(',')

    fig, axes = plt.subplots(
        nrows=2, ncols=2, figsize=(
            8 * 9, 4 * 9
        )
    )

    if axes.ndim == 1:
        axes = reshape(
            axes, (-1, 2)
        )

    fig.subplots_adjust(hspace=0.8)

    data = []
    columns = []
    additional = []
    seen = []

    if time:
        files = "*.time.data.tsv"
    else:
        files = "*.data.tsv"
    for file in directory.glob("*.data.tsv"):

        name = file.name.replace(files[1:], "")

        if nextflow:
            prefix_data = name.split(".")

            sketch = directory.name
            reads = prefix_data[-1]
            ranks = prefix_data[-2]
            name = ".".join(
                prefix_data[:-2]
            )

            additional.append(
                pandas.Series((sketch, ranks, reads))
            )

        data.append(
            pandas.read_csv(file, sep='\t', header=0)
        )
        columns.append(name)

    if data:
        all_data = pandas.concat(data, axis=1)
        features = all_data[['feature']]
        new_index = features.iloc[:, 0]

        predictions = all_data[['prediction']]
        predictions.index = new_index
        predictions.columns = columns

        stability = all_data[['stability']]
        stability.index = new_index
        stability.columns = columns

        preference = all_data[['preference']]
        preference.index = new_index
        preference.columns = columns

        if additional:
            additional_data = pandas.concat(additional, axis=1)
            additional_data.columns = columns
            additional_data.index = ['sketch', 'ranks', 'reads']

            predictions = pandas.concat([additional_data, predictions])
            stability = pandas.concat([additional_data, stability])
            preference = pandas.concat([additional_data, preference])

            if subset is not None:
                predictions, stability, preference = subset_data(
                    predictions, stability, preference, subset
                )

        if heatmap:

            for df in (predictions, stability, preference):
                df.drop(labels=['sketch', 'ranks', 'reads'], inplace=True)

            stability = stability.apply(pandas.to_numeric, errors='coerce')
            preference = preference.apply(pandas.to_numeric, errors='coerce')

            if threshold > 0:
                preference[preference < threshold] = 0

            plot_heatmap(
                values=preference,
                palette="Blues",
                ax=axes[0][0],
                threshold=0.,
                labels=predictions
            )

            plot_heatmap(
                values=preference,
                palette="Greens",
                ax=axes[0][1],
                threshold=0.6,
            )

            stability[stability == -1] = None

            plot_heatmap(
                values=stability+1,  # show reads not index
                palette="Purples_r",
                ax=axes[1][0],
                fmt=".0f",
            )

            plt.tight_layout()
            fig.savefig("test.png")

        predictions.T.sort_index().to_csv(f'{prefix}.predictions.tsv', sep='\t')
        stability.T.sort_index().to_csv(f'{prefix}.stability.tsv', sep='\t')
        preference.T.sort_index().to_csv(f'{prefix}.preference.tsv', sep='\t')
    else:
        print('No files found.')
        exit(1)


def plot_heatmap(
    values: pandas.DataFrame,
    palette: [str] or str = "YlGnBu",
    ax=None,
    fmt: str = ".2f",
    cbar: bool = True,
    annot: bool = True,
    labels: pandas.DataFrame = None,
    threshold: float = 0.
):

    # If all NA
    if all(values.isna().all().tolist()):
        values = values.fillna(0.)

    p1 = sns.heatmap(
        values.T.sort_index(),
        linewidths=5,
        cbar=cbar,
        ax=ax,
        annot=annot,
        fmt=fmt,
        cmap=palette,
        annot_kws={"size": 24, "weight": "bold"}
    )
    p1.tick_params(axis='both', which='major', labelsize=24, length=3, width=2)
    p1.tick_params(axis='x', rotation=90)
    p1.tick_params(axis='y', rotation=0)
    # p1.set_facecolor("lightgray")

    if threshold > 0.:
        for text in ax.texts:
            if float(text.get_text()) < threshold:
                text.set_text("")

    if labels is not None:
        label_vec = labels.T.sort_index().stack().tolist()
        for i, text in enumerate(ax.texts):
            try:
                # there is always categorical data never numeric floats
                val = f"{float(label_vec[i]):.0f}"
            except ValueError:
                val = label_vec[i].strip()
                if val.startswith('SCCmec'):
                    val = val.replace('SCCmec-', '')

            text.set_text(val)


def subset_data(predictions, stability, preference, subset: list):

    if len(subset) == 1:
        predictions = predictions.loc[
            predictions['']
        ]
    else:
        predictions = predictions.T
        predictions = predictions.loc[
            (predictions['ranks'] == subset[0]) &
            (predictions['reads'] == subset[1]), :
        ].T

        stability = stability.T
        stability = stability.loc[
            (stability['ranks'] == subset[0]) &
            (stability['reads'] == subset[1]), :
        ].T

        preference = preference.T
        preference = preference.loc[
            (preference['ranks'] == subset[0]) &
            (preference['reads'] == subset[1]), :
        ].T

    return predictions, stability, preference
