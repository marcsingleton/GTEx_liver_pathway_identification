"""Construct a network of differentially expressed liver genes."""

import os

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


cutoff = 0.75

df_corr = pd.read_table('../differential/out/pooled/corr.tsv', header=[0, 1], index_col=[0, 1])
gene_ids = df_corr.index.get_level_values('Name')
array = df_corr.to_numpy()

# Construct graph
graph = {gene_id: [] for gene_id in gene_ids}
for i in range(len(array)):
    gene_id1 = gene_ids[i]
    for j in range(i+1, len(array)):
        gene_id2 = gene_ids[j]
        r = array[i, j]
        if r >= cutoff:
            graph[gene_id1].append(gene_id2)
            graph[gene_id2].append(gene_id1)

# Identify connected components
components = get_connected_components(graph)

# Save results
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/graph.tsv', 'w') as file:
    file.write('node\tadjacent\n')
    for node, adjacents in graph.items():
        file.write(f'{node}\t{",".join(adjacents)}\n')

with open('out/components.tsv', 'w') as file:
    file.write('component_id\tgene_ids\n')
    for component_id, component in enumerate(components):
        file.write(f'{component_id:06}\t{",".join(component)}\n')
