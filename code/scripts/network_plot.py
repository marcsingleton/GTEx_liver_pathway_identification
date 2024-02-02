"""Plot statistics related to graph and clusters of differentially expressed liver genes."""

import os

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from matplotlib import cm
from matplotlib import colors
from matplotlib import ticker

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

# Plot n largest components
n = 10
fig_width = 4.8
margin_data = 0.025
cmap_base = plt.colormaps['plasma_r']
cmap = colors.ListedColormap(cmap_base.colors[int(0.25*len(cmap_base.colors)):])  # Trim lower 25% of colors
node_color = '#444444'
node_slope = 1000
node_intercept = 1
edge_slope = 10
edge_intercept = 0.75

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

    # Some formatting for node size and edges
    node_size = node_slope * xlen * ylen / len(graph_nx) + node_intercept  # Node size is inversely proportional to node density
    edge_width = edge_slope * xlen * ylen / len(graph_nx) + edge_intercept
    ws = [w for _, _, w in graph_nx.edges.data('weight')]
    wmin, wmax = min(ws), max(ws)

    # Make plot
    fig = plt.figure()
    ax = fig.add_axes((0, 0.075, 1, 0.925))
    nx.draw_networkx_edges(graph_nx, positions, ax=ax, width=edge_width, edge_color=ws, edge_cmap=cmap)
    nx.draw_networkx_nodes(graph_nx, positions, ax=ax, node_size=node_size, node_color=node_color, linewidths=0)

    ax.set_xlim((xmin - margin_data * xlen, xmax + margin_data * xlen))  # Set manually because draw_networkx_edges hard codes the data limits with 5% padding
    ax.set_ylim((ymin - margin_data * ylen, ymax + margin_data * ylen))
    ax.set_aspect('equal')
    ax.set_axis_off()

    cax = fig.add_axes((0.8, 0.05, 0.15, 0.015))
    cb = fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(wmin, wmax), cmap=cmap),
                      cax=cax, orientation='horizontal')
    cb.ax.xaxis.set_major_locator(ticker.LinearLocator(2))
    cb.ax.xaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.2}'))
    cb.outline.set_visible(False)

    fig.savefig(f'{prefix}/{rank:02}|{component_id}.png')
    plt.close()
