import click

import pandas
import dendropy
import networkx as nx

from pathlib import Path

from sketchy.sketchy import LineageIndex
from sketchy.utils import MutuallyExclusiveOption


@click.command()
@click.option(
    '--ssh',
    '-s',
    default=None,
    type=Path,
    required=True,
    help='Path to sum of shared hashes file to map to population [required]'
)
@click.option(
    '--tree',
    '-t',
    default=None,
    type=Path,
    required=False,
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["graph"],
    help='Path to phylogeny [newick] to use as population map [required]'
)
@click.option(
    '--graph',
    '-g',
    default=None,
    type=Path,
    required=False,
    cls=MutuallyExclusiveOption,
    mutually_exclusive=["tree"],
    help='Path to a population graph to use as population map [none]',
)
@click.option(
    '--index',
    '-i',
    default=None,
    type=Path,
    required=True,
    help='Path to db_lineage index file to use as validation population [required]'
)
@click.option(
    '--column',
    '-c',
    default='id',
    type=str,
    required=False,
    help='Column in index file that map the indices from the sum of '
         'shared hashes file to the population'
)
@click.option(
    '--output',
    '-o',
    default='pop.tsv',
    type=Path,
    required=False,
    help='Path to a output GIF [pop.gif]'
)
def popmap(ssh, tree, index, column, graph, output):

    """ Experimental: map reference sketch hits against population structures """

    lix = LineageIndex(index_file=index)

    if tree:
        t = dendropy.Tree.get(path=tree, schema="newick")
        genomes = [
            str(tax).replace("'", "") for tax in t.taxon_namespace
        ]
    elif graph:
        g = nx.read_edgelist(str(graph))
        genomes = [
            str(node) for node, data in g.nodes(data=True)
        ]
        g = nx.relabel.convert_node_labels_to_integers(g)
        nx.write_edgelist(g, str(graph.stem)+'.iel', data=False)
    else:
        raise ValueError(
            'Either a population graph (--graph) or phylogenetic tree (--tree) '
            'are required for population mapping.'
        )

    ssh_reads = pandas.read_csv(
        ssh, sep='\t', names=['idx', 'ssh', 'rank', 'read']
    )

    resolved = ssh_reads.merge(lix.index, how='left', on='idx')

    genomes_by_read = [
        [1 if g in read_data[column].values else 0 for g in genomes]
        for i, read_data in resolved.groupby('read', sort=False)
    ]

    df = pandas.DataFrame(data=genomes_by_read)
    df.columns = genomes
    df.to_csv(output, sep='\t', index=False)