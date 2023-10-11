"""Initial exploratory analyses of GTEx data for identification of liver-specific pathways.

The main goal of this script is to perform an exploratory analysis of the GTEx
data, hopefully identifying any outliers or other unusual aspects of the data
that may need special treatment before a subsequent differential expression
analysis.
"""

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from src.constants import tissue_regex

data_path_read = '../../data/GTEx/bulk_count/'
data_path_gene = '../protein_coding/out/genes.tsv'

# Load gene metadata
df_gene = pd.read_table(data_path_gene)
genes_protein = df_gene.loc[df_gene['gene_type'] == 'protein_coding', 'gene_id']
genes_other = df_gene.loc[df_gene['gene_type'] != 'protein_coding', 'gene_id']

# Load mapped reads
dfs_read = []
files = sorted([file for file in os.listdir(data_path_read) if file.endswith('.gct.gz')])
for file in files:
    # Read in data and transpose to make genes columns
    df_read = pd.read_table(f'{data_path_read}/{file}', skiprows=2)
    df_read = df_read.drop('id', axis=1).set_index(['Name', 'Description'])
    df_read.columns.name = 'Sample'
    df_read = df_read.transpose()

    # Add tissue information to index
    tissue = re.search(tissue_regex, file).group(1)
    df_read['Tissue'] = tissue
    df_read = df_read.set_index('Tissue', append=True)

    # Drop non-protein coding genes
    df_read = df_read.drop(genes_other, axis=1, level=0)
    if set(df_read.columns.get_level_values('Name')) != set(genes_protein):
        print('Warning: gene_ids in read dataframe do not match protein-coding gene_ids.')

    dfs_read.append(df_read)
df_read = pd.concat(dfs_read)

prefix = 'out/'
if not os.path.exists(prefix):
    os.mkdir(prefix)

# Bar graph of number of samples for each tissue type
counts = df_read.groupby(level='Tissue').size()
xs = list(range(len(counts)))
ys = counts.values
labels = counts.index

fig, ax = plt.subplots(gridspec_kw={'bottom': 0.25})
ax.bar(xs, ys)
ax.set_xticks(xs, labels, rotation=45, horizontalalignment='right', fontsize=8)
ax.set_ylabel('Number of samples')
fig.savefig(f'{prefix}/bar|sample_number-tissue.png')
plt.close()

# Bar graph of mean number of reads by tissue
mean = df_read.sum(axis=1).groupby(level='Tissue').mean()
xs = list(range(len(mean)))
ys = mean.values
labels = mean.index

fig, ax = plt.subplots(gridspec_kw={'bottom': 0.25})
ax.bar(xs, ys)
ax.set_xticks(xs, labels, rotation=45, horizontalalignment='right', fontsize=8)
ax.set_ylabel('Mean reads in samples')
fig.savefig(f'{prefix}/bar|reads_mean-tissue.png')
plt.close()

# Bar graph of number of unique donors by tissue
index = df_read.index.to_frame(index=False)
index['Donor'] = index['Sample'].apply(lambda x: x.split('-')[1])
counts = index.groupby('Tissue')['Donor'].nunique()
xs = list(range(len(counts)))
ys = counts.values
labels = counts.index

fig, ax = plt.subplots(gridspec_kw={'bottom': 0.25})
ax.bar(xs, ys)
ax.set_xticks(xs, labels, rotation=45, horizontalalignment='right', fontsize=8)
ax.set_ylabel('Number of unique donors')
fig.savefig(f'{prefix}/bar|donor_number-tissue.png')
plt.close()

# Violin plot of reads by tissue
groups = [(label, df.to_numpy()) for label, df in df_read.sum(axis=1).groupby(level='Tissue')]
labels, arrays = list(zip(*groups))
positions = list(range(len(groups)))

fig, ax = plt.subplots(gridspec_kw={'bottom': 0.25})
ax.violinplot(arrays, positions=positions, showmedians=True)
ax.set_xticks(positions, labels, rotation=45, horizontalalignment='right', fontsize=8)
ax.set_ylabel('Number of reads in sample')
fig.savefig(f'{prefix}/violin|reads-tissue.png')
plt.close()

# Histogram plot of gene count z-scores by tissue
eps = 1E-6  # Use a small epsilon for the denominator to prevent division by zero for constant genes
zscore = df_read.groupby(level='Tissue').transform(lambda x: (x - x.mean()) / (x.std() + eps))
groups = [(label, df.to_numpy().ravel()) for label, df in zscore.groupby(level='Tissue')]
labels, arrays = list(zip(*groups))

fig, axs = plt.subplots(2, len(groups) // 2, figsize=(12.8, 4.8), layout='constrained')
for ax, label, array in zip(axs.ravel(), labels, arrays):
    ax.hist(array, bins=100)
    ax.tick_params(labelsize=6)
    ax.ticklabel_format(scilimits=(-6, 6))
    offset_text = ax.yaxis.get_offset_text()
    offset_text.set_fontsize(6)
    offset_text.set_position((-0.1, 1))
    ax.set_title(label, fontsize=8)
fig.supxlabel('$z$-score of reads to gene in tissue and gene group', fontsize=10)
fig.supylabel('Number of genes', fontsize=10)
fig.savefig(f'{prefix}/hist|gene_number-zscore.png')
plt.close()

# PCA plots of samples with gene counts as features
pca = PCA(n_components=10)
transform = pca.fit_transform(df_read.to_numpy())
labels = df_read.index.get_level_values('Tissue')

# Scatters of projection onto PCs
fig, ax = plt.subplots(gridspec_kw={'right': 0.7})
for label in labels.unique():
    ax.scatter(transform[label == labels, 0], transform[label == labels, 1],
               s=8, alpha=0.75, edgecolor='none', label=label)
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8, markerscale=2)
fig.savefig(f'{prefix}/scatter|PC1-PC2.png')
plt.close()

fig, ax = plt.subplots(gridspec_kw={'right': 0.7})
for label in labels.unique():
    ax.scatter(transform[label == labels, 1], transform[label == labels, 2],
               s=8, alpha=0.75, edgecolor='none', label=label)
ax.set_xlabel('PC2')
ax.set_ylabel('PC3')
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8, markerscale=2)
fig.savefig(f'{prefix}/scatter|PC2-PC3.png')
plt.close()

# Bar graph of explained variance ratio
ys = pca.explained_variance_ratio_
xs = list(range(1, len(ys)+1))

fig, ax = plt.subplots()
ax.bar(xs, ys)
ax.set_xlabel('PC')
ax.set_ylabel('Explained variance ratio')
fig.savefig(f'{prefix}/bar|scree.png')
plt.close()
