import click
import pandas

from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns


@click.command()
@click.option(
    '--report', '-r', type=str, help='Path or file glob to tax report files'
)
@click.option(
    '--prefix', '-p', type=str, help='Output prefix for plot file.'
)
@click.option(
    '--level', '-l', type=str, default='S',
    help='Taxonomic level to assess: species [S]'
)
@click.option(
    '--top', '-t', type=int, default=10,
    help='Show top taxonomic levels in plots [10]'
)
@click.option(
    '--color', '-c', type=str, default='Greens',
    help='Color palette for central donut plot.'
)
@click.option(
    '--title', '-t', type=str, default=None,
    help='Row titles for center plot, comma separated string.'
)
@click.option(
    '--sub', '-s', is_flag=True,
    help='Add subplot titles for each column.'
)
def kraken_report(report, prefix, level, top, color, title, sub):
    """ Generate metagenomic report plots from Kraken2 in Sketchy Nextflow """

    if "*" in report:
        reports = list(Path().glob(
            str(report)
        ))
    elif "," in report:
        reports = [Path(p) for p in report.split(',')]
    else:
        reports = [Path(report)]

    fig, axes = plt.subplots(
        nrows=len(reports), ncols=3, figsize=(27.0, len(reports) * 9),
        subplot_kw=dict(aspect="equal")
    )

    for i, report in enumerate(reports):

        df = pandas.read_csv(
            report, header=None, sep="\t", names=[
                "percent", "reads", "direct", "level", "taxid", "taxonomy",
            ]
        )
        df.taxonomy = df.taxonomy.str.strip()

        if len(reports) == 1:
            (ax1, ax2, ax3) = axes
        else:
            (ax1, ax2, ax3) = axes[i, :]
        # Plots

        human, unclassified, ureads = plot_overview(df=df, ax=ax1, color=color)

        creads = int(
            df[df.taxonomy == 'root'].reads
        )

        total = int(creads + ureads)

        plot_major_composition(
            df=df, ax=ax2, level=level, other=human+unclassified, color=color
        )
        plot_minor_composition(df=df, ax=ax3, level=level, top=top)

        if sub:
            ax1.title.set_text(f'Reads ({total})')
            ax2.title.set_text(f'Major Taxa (> 10%)')
            ax3.title.set_text('Minor Taxa (<= 10%)')
        # Final output

        plt.tight_layout()
        fig.subplots_adjust(top=1)

        fig.suptitle(title, fontsize=24)

        fig.savefig(f'{prefix}.svg', pad_inches=0.5)
        fig.savefig(f'{prefix}.pdf', pad_inches=0.5)

def plot_overview(df, ax, color=None):

    target = sns.color_palette(color, 2)[-1]

    overview_data, overview_labels = [], []

    human = df[df['taxonomy'] == 'Homo sapiens']

    if not human.empty:
        overview_data.append(
            float(human.percent)
        )
        overview_labels.append(f'Human [{float(human.percent)}%]')

    unclassified = df[df['taxonomy'] == 'unclassified']

    if not unclassified.empty:
        overview_data.append(
            float(unclassified.percent)
        )
        overview_labels.append(f'Unclassified [{float(unclassified.percent)}%]')

    other = round(100 - sum(overview_data), ndigits=2)

    overview_data.append(other)
    overview_labels.append(f'Targets [{other}%]')

    if len(overview_labels) == 2:
        # No human
        colors = ['#A9A9A9', target]
    else:
        # With human
        colors = ['#fec44f', '#A9A9A9', target]

    plot_annotated_donut(
        labels=overview_labels, data=overview_data, ax=ax, colors=colors
    )

    if human.empty:
        human = 0
    else:
        human = human.percent

    if unclassified.empty:
        un = 0
    else:
        un = unclassified.percent

    return float(human), float(un), int(unclassified.reads)


def plot_major_composition(df, ax, level, other=0., color=None):

    df2 = df[(~df.taxonomy.isin(
        ['Homo sapiens', 'Homo']
    )) & (df.percent >= 10.)]

    tax = df2[
        df2['level'] == level
    ].sort_values(by='percent', ascending=False)

    data = tax.percent.tolist()
    labels = tax.taxonomy.tolist()

    minor = round(100 - other - sum(data), ndigits=2)

    data.insert(0, minor)
    labels.insert(0, 'Minor')

    labels = [f'{l}: {data[i]}%' for i, l in enumerate(labels)]

    color = sns.color_palette(color, len(data))
    color[0] = '#808080'

    plot_annotated_donut(
        labels=labels, data=data, ax=ax, colors=color
    )


def plot_minor_composition(df, ax, level, top):

    # Remove human and filter percent < 10%

    df2 = df[(~df.taxonomy.isin(
        ['Homo sapiens', 'Homo']
    )) & (df.percent < 10.)]

    tax = df2[
        df2['level'] == level
    ].sort_values(by='percent', ascending=False)

    df = tax[:top]

    df = df.sort_values(by='percent', ascending=True)

    df.plot(
        y='percent', x='taxonomy', kind='barh', ax=ax, legend=False, color='#808080'
    )

    values = df.percent.tolist()
    ax.set_xlabel('Percent of reads classified (%)')
    ax.set_ylabel('')
    ax.set_xlim(0, 10)

    # set individual bar lables using above list
    for i, patch in enumerate(ax.patches):
        # get_width pulls left or right; get_y pushes up or down
        ax.text(
            patch.get_width() + .3, patch.get_y() + 0.1,
            str(values[i]) + "%",
            fontsize=10,
            color='dimgrey'
        )


def plot_annotated_donut(labels: [str], data: [float], ax: None, colors=None):

    """ Donut chart with labels example
    """

    if len(labels) != len(data):
        raise ValueError('Data and label lists most be of the same length.')

    wedges, texts = ax.pie(
        data, wedgeprops=dict(width=0.5), startangle=-40, colors=colors
    )

    bbox_props = dict(
        boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72
    )

    kw = dict(
        arrowprops=dict(arrowstyle="-"),
        bbox=bbox_props, zorder=0, va="center")

    for i, p in enumerate(wedges):
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = "angle,angleA=0,angleB={}".format(ang)
        kw["arrowprops"].update({"connectionstyle": connectionstyle})
        ax.annotate(labels[i], xy=(x, y), xytext=(1.35 * np.sign(x), 1.4 * y),
                    horizontalalignment=horizontalalignment, **kw, fontsize=16)