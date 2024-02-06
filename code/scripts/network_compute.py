"""Construct a network of differentially expressed liver genes."""

import argparse
import os

import numpy as np
import pandas as pd


def get_connected_components(graph):
    """Find connected components of graph, implemented as a depth-first search."""
    # Initialize
    component = None
    components = []
    expand_stack = []  # Stack to expand current component
    search_stack = sorted(graph)  # Stack to search for new component
    marked = set()

    # Loop
    while expand_stack or search_stack or component is not None:
        # Exhaust expand stack first
        while expand_stack:
            node = expand_stack.pop()
            if node in marked:
                continue
            component.add(node)
            expand_stack.extend(sorted(graph[node]))
            marked.add(node)
        if component is not None:  # Only record component if not None; only False in first iteration
            components.append(component)
            component = None

        # Proceed to search stack when expand stack is empty
        while search_stack and component is None:
            node = search_stack.pop()
            if node in marked:  # Skip previously added nodes; necessary to prevent loops and incomplete components
                continue
            component = {node}
            expand_stack.extend(sorted(graph[node]))
    return components

parser = argparse.ArgumentParser()
parser.add_argument('corr')
parser.add_argument('-o', '--output_path', default='./')

top_N = 4

if __name__ == '__main__':
    args = parser.parse_args()

    df_corr = pd.read_table(args.corr, header=[0, 1], index_col=[0, 1])
    gene_ids = df_corr.index.get_level_values('Name')
    array = df_corr.to_numpy()

    # Make sets of top correlations for each gene
    inf_diag = np.copy(array)
    np.fill_diagonal(inf_diag, -np.inf)
    array_sorted = -np.sort(-inf_diag, axis=1)
    index_sorted = np.argsort(-inf_diag, axis=1)

    top_genes_map = {}
    for i in range(len(array_sorted)):
        r_cutoff = array_sorted[i, top_N-1]
        gene_id1 = gene_ids[i]
        top_genes = set()
        for j in range(len(array_sorted)):
            r = array_sorted[i, j]
            if r >= r_cutoff:
                index = index_sorted[i, j]
                gene_id2 = gene_ids[index]
                top_genes.add(gene_id2)
            else:
                break
        top_genes_map[gene_id1] = top_genes

    # Construct graph
    graph = {gene_id: [] for gene_id in gene_ids}
    for i in range(len(array)):
        gene_id1 = gene_ids[i]
        top_genes1 = top_genes_map[gene_id1]
        for j in range(i+1, len(array)):
            gene_id2 = gene_ids[j]
            top_genes2 = top_genes_map[gene_id2]
            if gene_id1 in top_genes2 and gene_id2 in top_genes1:
                graph[gene_id1].append(gene_id2)
                graph[gene_id2].append(gene_id1)

    # Identify connected components
    components = get_connected_components(graph)

    # Save results
    prefix = args.output_path
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    with open(f'{prefix}/graph.tsv', 'w') as file:
        file.write('node\tadjacent\n')
        for node, adjacents in graph.items():
            file.write(f'{node}\t{",".join(adjacents)}\n')

    with open(f'{prefix}/components.tsv', 'w') as file:
        file.write('component_id\tgene_ids\n')
        for component_id, component in enumerate(components):
            file.write(f'{component_id:06}\t{",".join(component)}\n')
