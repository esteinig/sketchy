import click
from numpy import where

from pathlib import Path
from sketchy.sketchy import SketchyDiagnostics


@click.command()
@click.option(
    '--nextflow',
    '-n',
    type=Path,
    required=True,
    help='Path to collected results directory of Nextflow run [required]',
)
@click.option(
    '--outdir',
    '-o',
    type=Path,
    default=Path("nfx-heatmaps"),
    help='Output directory for plots [nfx-heatmaps]'
)
@click.option(
    '--plot',
    '-p',
    type=Path,
    default=Path("diagnostics.png"),
    help='Plot file, extension specifies format [diagnostics.png]'
)
@click.option(
    '--color',
    '-c',
    type=str,
    default='YlGnBu',
    help='Color palette for output plots [YlGnBu]'
)
@click.option(
    '--subset_column',
    '-sc',
    type=str,
    default=None,
    help='Column name to subset data [none]'
)
@click.option(
    '--subset_values',
    '-sv',
    type=str,
    default='None',
    help='Comma delimited string of values in subset column to select [none]'
)
@click.option(
    '--verbose',
    '-v',
    is_flag=True,
    help='Verbose logging output [false]'
)
@click.option(
    '--mpl_backend',
    type=str,
    default="",
    help='Matplotlib backend [default]'
)
def genotype_heatmap(nextflow, outdir, plot, color, subset_column, subset_values, mpl_backend, verbose):

    """ Comparison of genotype predictions from Nextflow """

    sd = SketchyDiagnostics(outdir=outdir, verbose=verbose, mpl_backend=mpl_backend)

    sd.plot_genotype_heatmap(nextflow=nextflow, subset_column=subset_column, subset_values=subset_values)

#
#     if subset:
#         subset = subset.split(',')
#
#     if time:
#         nrows = 2
#         ncols = 2
#     else:
#         nrows = 3
#         ncols = 1
#
#     if reference is not None and time:
#         nrows = 3
#         ncols = 2
#     elif reference is not None and not time:
#         nrows = 5
#         ncols = 1
#
#     fig, axes = plt.subplots(
#         nrows=nrows, ncols=ncols, figsize=(
#             ncols * 4 * 9 * scale, nrows * 2 * 9 * scale
#         )
#     )
#
#     data = []
#     columns = []
#     additional = []
#
#     if time:
#         files = "*.time.data.tsv"
#     else:
#         files = "*.data.tsv"
#
#     if statistics:
#         files = "*.filtered.stats.txt"
#
#     if coverage:
#         files = "*.coverage.txt"
#
#     stats = []
#     covs = []
#     for file in directory.glob(files):
#         if statistics:
#             s = pandas.read_csv(
#                 file,
#                 sep=' ',
#                 names=[
#                     'reads',
#                     'bp',
#                     'max_length',
#                     'min_length',
#                     'mean_length',
#                     'median_length',
#                     'mean_quality',
#                     'median_quality'
#                 ]
#             )
#
#             s['prefix'] = file.name.replace(files[1:], '')
#             stats.append(s)
#             continue
#
#         if time and 'time' not in file.name:
#             continue
#         if not time and 'time' in file.name:
#             continue
#
#         name = file.name.replace(files[1:], "")
#
#         if nextflow:
#             prefix_data = name.split(".")
#
#             sketch = directory.name
#             reads = prefix_data[-1]
#             ranks = prefix_data[-2]
#             name = ".".join(
#                 prefix_data[:-2]
#             )
#
#             additional.append(
#                 pandas.Series((sketch, ranks, reads))
#             )
#
#         data.append(
#             pandas.read_csv(file, sep='\t', header=0)
#         )
#         columns.append(name)
#
#     if data:
#         predictions, stability, preference, times = get_data(data, columns)
#
#         if additional:
#             additional_data = pandas.concat(additional, axis=1)
#             additional_data.columns = columns
#             additional_data.index = ['sketch', 'ranks', 'reads']
#
#             predictions = pandas.concat([additional_data, predictions])
#             stability = pandas.concat([additional_data, stability])
#             preference = pandas.concat([additional_data, preference])
#             times = pandas.concat([additional_data, times])
#
#             if subset is not None:
#                 predictions, stability, preference, times = subset_data(
#                     predictions, stability, preference, times, subset
#                 )
#
#         predictions = predictions.T.sort_index()
#         stability = stability.T.sort_index()
#         preference = preference.T.sort_index()
#
#         if times is not None:
#             times = times.T.sort_index()
#
#         if heatmap:
#
#             for df in (predictions, stability, preference, times):
#                 try:
#                     df.drop(columns=['sketch', 'ranks', 'reads'], inplace=True)
#                 except (KeyError, AttributeError):
#                     continue
#
#             stability = stability.apply(pandas.to_numeric, errors='coerce')
#             preference = preference.apply(pandas.to_numeric, errors='coerce')
#
#             if threshold > 0:
#                 preference[preference < threshold] = 0
#
#             plot_heatmap(
#                 values=preference,
#                 palette="Greens",
#                 ax=axes[0][0] if time else axes[0],
#                 threshold=0.,
#                 labels=predictions,
#                 title="\nSketchy predictions\n"
#             )
#
#             plot_heatmap(
#                 values=preference,
#                 palette="BuGn",
#                 ax=axes[0][1] if time else axes[1],
#                 threshold=threshold,
#                 title="\nMedian preference score\n"
#             )
#
#             stability[stability == -1] = None
#
#             plot_heatmap(
#                 values=stability,  # show reads not index
#                 palette="Blues_r",
#                 ax=axes[1][0] if time else axes[2],
#                 fmt=".0f",
#                 title="\nStable breakpoint (reads)\n"
#             )
#
#             if time:
#                 plot_heatmap(
#                     values=stability,  # show reads not index
#                     palette="PuBu_r",
#                     ax=axes[1][1],
#                     fmt=".0f",
#                     labels=times,
#                     time=True,
#                     title="\nStable breakpoint (time)\n"
#                 )
#
#             if reference is not None:
#                 reference_genotypes = pandas.read_csv(
#                     reference, sep='\t', header=0, index_col=0
#                 )
#
#                 evaluation = compare_to_prediction(
#                     reference_genotypes, predictions, preference, threshold
#                 )
#
#                 plot_heatmap(
#                     values=evaluation,
#                     palette="PiYG_r",
#                     ax=axes[2][0] if time else axes[3],
#                     threshold=0.,
#                     labels=reference_genotypes,
#                     evaluation=True,
#                     title="\nReference + evaluation\n"
#                 )
#
#                 plot_heatmap(
#                     values=evaluation,
#                     palette="PRGn_r",
#                     ax=axes[2][1] if time else axes[4],
#                     threshold=0.,
#                     labels=predictions,
#                     evaluation=True,
#                     title="\nPrediction + evaluation\n"
#                 )
#
#             plt.tight_layout()
#             fig.savefig(f"{prefix}.{image_format}")
#
#         predictions.to_csv(f'{prefix}.predictions.tsv', sep='\t')
#         stability.to_csv(f'{prefix}.stability.tsv', sep='\t')
#         preference.to_csv(f'{prefix}.preference.tsv', sep='\t')
#     else:
#         print('No files found.')
#         exit(1)
#
#
# def plot_heatmap(
#         values: pandas.DataFrame,
#         palette: [str] or str = "YlGnBu",
#         ax=None,
#         fmt: str = ".3f",
#         cbar: bool = True,
#         annot: bool = True,
#         labels: pandas.DataFrame = None,
#         threshold: float = 0.,
#         time: bool = False,
#         title: str = "",
#         evaluation: bool = False
# ):
#     # If all NA
#     if all(
#             values.isna().all().tolist()
#     ):
#         values = values.fillna(0.)
#
#     p1 = sns.heatmap(
#         values,
#         vmin=0 if evaluation else None,
#         vmax=3 if evaluation else None,
#         linewidths=5,
#         cbar=cbar,
#         ax=ax,
#         annot=annot,
#         fmt=fmt,
#         cmap=palette,
#         annot_kws={"size": 18 if time else 24, "weight": "bold"}
#     )
#     p1.tick_params(axis='both', which='major', labelsize=24, length=3, width=2)
#     p1.tick_params(axis='x', rotation=90)
#     p1.tick_params(axis='y', rotation=0)
#
#     ax.set_title(title, fontdict={'fontsize': 24})
#     ax.set_xlabel('')
#
#     if threshold > 0.:
#         for text in ax.texts:
#             if float(text.get_text()) < threshold:
#                 text.set_text("")
#
#     if not time and labels is not None:
#         label_vec = labels.stack().tolist()
#         for i, text in enumerate(ax.texts):
#             try:
#                 # there is always categorical data never numeric floats
#                 val = f"{float(label_vec[i]):.0f}"
#             except ValueError:
#                 val = label_vec[i].strip()
#                 if val.startswith('SCCmec'):
#                     val = val.replace('SCCmec-', '')
#
#             text.set_text(val)
#
#     if time:
#         label_vec = labels.stack().tolist()
#         for i, text in enumerate(ax.texts):
#             if isinstance(label_vec[i], str):
#                 try:
#                     text.set_text(
#                         label_vec[i].split(" ")[1]
#                     )  # date time
#                 except KeyError:
#                     text.set_text(
#                         label_vec[i]
#                     )  # time delta
#             else:
#                 text.set_text()
#
#
# def get_data(data: list, columns: list):
#     all_data = pandas.concat(data, axis=1)
#
#     features = all_data[['feature']]
#     new_index = features.iloc[:, 0]
#
#     predictions = all_data[['prediction']]
#     predictions.index = new_index
#     predictions.columns = columns
#
#     stability = all_data[['stability']]
#     stability.index = new_index
#     stability.columns = columns
#
#     preference = all_data[['preference']]
#     preference.index = new_index
#     preference.columns = columns
#
#     if 'time' in all_data.columns:
#         time = all_data[['time']]
#         time.index = new_index
#         time.columns = columns
#     else:
#         time = None
#
#     return predictions, stability, preference, time
#
#
# def subset_data(predictions, stability, preference, times, subset: list):
#     if len(subset) == 1:
#         predictions = predictions.loc[
#             predictions['']
#         ]
#     else:
#         predictions = predictions.T
#         predictions = predictions.loc[
#                       (predictions['ranks'] == subset[0]) &
#                       (predictions['reads'] == subset[1]), :
#                       ].T
#
#         stability = stability.T
#         stability = stability.loc[
#                     (stability['ranks'] == subset[0]) &
#                     (stability['reads'] == subset[1]), :
#                     ].T
#
#         preference = preference.T
#         preference = preference.loc[
#                      (preference['ranks'] == subset[0]) &
#                      (preference['reads'] == subset[1]), :
#                      ].T
#
#         times = times.T
#         times = times.loc[
#                 (times['ranks'] == subset[0]) &
#                 (times['reads'] == subset[1]), :
#                 ].T
#
#     return predictions, stability, preference, times
#
#
# def compare_to_prediction(reference, prediction, preference, threshold):
#     """ Compare reference to predictions and evaluate
#
#     :param reference:
#     :param prediction:
#     :return:
#
#     """
#
#     positives = reference == prediction
#     negatives = reference != prediction
#
#     evaluation = reference.copy()
#
#     confident_positives = positives & (preference >= threshold)
#     unconfident_positives = positives & (preference < threshold)
#
#     confident_negatives = negatives & (preference >= threshold)
#     unconfident_negatives = negatives & (preference < threshold)
#
#     # from good to bad essentially
#     for i, mask in enumerate((
#             confident_positives,
#             unconfident_positives,
#             unconfident_negatives,
#             confident_negatives,
#     )):
#         evaluation[mask] = i
#
#     return evaluation.apply(pandas.to_numeric, errors='coerce')