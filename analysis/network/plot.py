"""Plot statistics related to graph and clusters of differentially expressed liver genes."""

import os

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

data_path_corr = '../differential/out/pooled/corr.tsv'
data_path_graph = 'out/graph.tsv'
data_path_components = 'out/components.tsv'

# Load graph
graph = {}
with open(data_path_graph) as file:
    file.readline()  # Skip header
    for line in file:
        node, adjacents = line.rstrip('\n').split('\t')
        adjacents = adjacents.split(',')
        graph[node] = adjacents

# Load components
components = {}
with open(data_path_components) as file:
    file.readline()  # Skip header
    for line in file:
        component_id, gene_ids = line.rstrip('\n').split('\t')
        gene_ids = gene_ids.split(',')
        components[component_id] = gene_ids

# Calculate statistics related to components
records = []
for component_id, gene_ids in components.items():
    records.append({'component_id': component_id,
                    'gene_num': len(gene_ids),
                    'edge_num': sum([len(graph[gene_id]) for gene_id in gene_ids])})
df = pd.DataFrame(records)

prefix = 'out/'

# Bar graph of number of genes in components
counts = df['gene_num'].value_counts()

fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of genes in component')
ax.set_ylabel('Number of components')
fig.savefig(f'{prefix}/bar|component_number-gene_number.png')
plt.close()

# Plot largest n components
n = 10
margin_data = 0.01

prefix = 'out/networks/'
if not os.path.exists(prefix):
    os.makedirs(prefix)

component_ids = df.sort_values(by='gene_num', ascending=False, ignore_index=True).loc[:n, 'component_id']
for component_id in component_ids:
    component = components[component_id]

    # Construct graph
    graph_nx = nx.Graph()
    for node in component:
        graph_nx.add_node(node)
        for adjacent in graph[node]:
            if (node, adjacent) not in graph_nx.edges:
                graph_nx.add_edge(node, adjacent)

    # Get positions and axes limits
    positions = nx.kamada_kawai_layout(graph_nx)
    xs = [xy[0] for xy in positions.values()]
    xmin, xmax = min(xs), max(xs)
    xlen = xmax - xmin
    ys = [xy[1] for xy in positions.values()]
    ymin, ymax = min(ys), max(ys)
    ylen = ymax - ymin

    # Make plot
    fig, ax = plt.subplots(layout='constrained')
    nx.draw_networkx_edges(graph_nx, positions, ax=ax, width=0.5, edge_color='black', alpha=0.1)
    nx.draw_networkx_nodes(graph_nx, positions, ax=ax, node_size=4, node_color='C0', linewidths=0)
    ax.set_xlim((xmin - margin_data * xlen, xmax + margin_data * xlen))  # Set manually because draw_networkx_edges hard codes the data limits with 5% padding
    ax.set_ylim((ymin - margin_data * ylen, ymax + margin_data * ylen))
    ax.set_aspect('equal')
    ax.set_axis_off()
    fig.savefig(f'{prefix}/{component_id}.png')
    plt.close()
