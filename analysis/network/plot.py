"""Plot statistics related to graph and clusters of differentially expressed liver genes."""

import os

import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pandas as pd
from matplotlib import colors
from matplotlib import cm

data_path_corr = '../differential/out/pooled/corr.tsv'
data_path_graph = 'out/graph.tsv'
data_path_components = 'out/components.tsv'

df_corr = pd.read_table(data_path_corr, header=[0, 1], index_col=[0, 1])
gene_ids = df_corr.index.get_level_values('Name')
indexes = {gene_id: index for index, gene_id in enumerate(gene_ids)}
array = df_corr.to_numpy()

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
ax.set_yscale('log')
fig.savefig(f'{prefix}/bar|component_number-gene_number|log.png')
plt.close()

# Plot largest n components
n = 10
fig_width = 4.8
margin_data = 0.05
cmap_base = plt.colormaps['plasma_r']
cmap = colors.ListedColormap(cmap_base(np.linspace(0.25, 1, 256)))
node_color = '#444444'

prefix = 'out/networks/'
if not os.path.exists(prefix):
    os.makedirs(prefix)

component_ids = df.sort_values(by='gene_num', ascending=False, ignore_index=True).loc[:n, 'component_id']
for rank, component_id in enumerate(component_ids):
    component = components[component_id]

    # Construct graph
    graph_nx = nx.Graph()
    for node1 in component:
        graph_nx.add_node(node1)
        for node2 in graph[node1]:
            if (node1, node2) not in graph_nx.edges:
                index1 = indexes[node1]
                index2 = indexes[node2]
                graph_nx.add_edge(node1, node2, weight=array[index1, index2])

    # Get positions and axes limits
    positions = nx.kamada_kawai_layout(graph_nx)
    xs = [xy[0] for xy in positions.values()]
    xmin, xmax = min(xs), max(xs)
    xlen = xmax - xmin
    ys = [xy[1] for xy in positions.values()]
    ymin, ymax = min(ys), max(ys)
    ylen = ymax - ymin

    # Rotate positions
    if xlen / ylen < 1:  # Make width longer side
        xmin, xmax, ymin, ymax = ymin, ymax, -xmax, -xmin
        xlen, ylen = ylen, xlen
        positions = {node: (y, -x) for node, (x, y) in positions.items()}

    # Adjust dimensions so aspect ratio is 1:1
    fig_height = fig_width * ylen / xlen
    figsize = (fig_width, fig_height)

    # Draw graph labeled by edge
    edges = sorted(graph_nx.edges, key=lambda x: graph_nx.edges[x]['weight'])
    ws = [w for _, _, w in graph_nx.edges.data('weight')]
    wmin, wmax = min(ws), max(ws)
    wlen = wmax - wmin

    # Make plot
    fig, ax = plt.subplots(layout='constrained')
    nx.draw_networkx_edges(graph_nx, positions, ax=ax, width=0.75, edge_color=ws, edge_cmap=cmap)
    nx.draw_networkx_nodes(graph_nx, positions, ax=ax, node_size=12, node_color=node_color, linewidths=0)

    ax.set_xlim((xmin - margin_data * xlen, xmax + margin_data * xlen))  # Set manually because draw_networkx_edges hard codes the data limits with 5% padding
    ax.set_ylim((ymin - margin_data * ylen, ymax + margin_data * ylen))
    ax.set_axis_off()

    ticks = [wmin + wlen / 4, wmax - wlen / 4]
    ticklabels = [f'{tick:.2}' for tick in ticks]
    cax = ax.inset_axes((0.8, 0.01, 0.15, 0.015))
    cbar = fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(wmin, wmax), cmap=cmap), cax=cax, orientation='horizontal')
    cbar.ax.set_xticks(ticks, ticklabels, fontsize=6)
    cbar.outline.set_visible(False)

    fig.savefig(f'{prefix}/{rank:02}|{component_id}.png')
    plt.close()
