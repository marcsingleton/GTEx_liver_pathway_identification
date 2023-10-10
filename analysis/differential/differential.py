"""Identification of differentially expressed genes in the liver.

The overall approach for identifying differentially expressed genes in the
liver was originally a simple one-way ANOVA on the tissue type where each gene
is treated separately. Because under this analysis, essentially all genes were
differentially expressed, I instead pivoted to a more straightforward
comparison of two means approach where the liver samples formed the
experimental group and all other samples were pooled together as the background
or "control" group. I assumed the variances were unequal based on the biology
of the two groups, and I also used a one-sided "greater" alternative to select
up-regulated genes only.

More sophisticated analyses can use more complex tools like mixed
effects models or DESeq2, which pool information across samples and/or genes to
obtain better estimates of key parameters.
"""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import f_oneway, ttest_ind


def benjamini_hochberg(pvalues, q):
    """Return array of BH procedure significance results.

    True indicates test is significant at the given FDR (q). Results are sorted
    by order of p-values in the input array."""
    pvalues = np.array(pvalues)

    # Sort values and track original indices
    ranked = np.sort(pvalues)
    index = np.argsort(pvalues)

    # Remove invalid tests and create ranks
    b = ~np.isnan(ranked)
    ranked = ranked[b]
    index = index[b]
    ranks = np.arange(1, len(ranked)+1)

    # BH procedure
    bh_critical = ranks / len(ranks) * q
    bh_test = ranked <= bh_critical
    bh_index = len(ranks) - 1 - np.argmax(np.flip(bh_test))  # argmax finds first true value; flip array to find last

    # Report results
    result = np.full(len(pvalues), False)
    result[index[:bh_index+1]] = True  # Flip all indices in original array corresponding to those below the cutoff to True

    return result


data_path_read = '../../data/GTEx/bulk_count/'
data_path_gene = '../protein_coding/out/genes.tsv'
tissue_regex = r'v8_([a-z_]+)\.gct\.gz'

q_BH = 0.01  # BH false discovery rate

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

# Normalize each tissue by read depth
total_reads = df_read.sum(axis=1)
df_read = df_read.mul(1 / total_reads, axis=0)

prefix = 'out/unpooled/'
if not os.path.exists(prefix):
    os.makedirs(prefix)

# Calculate ANOVAs (unpooled)
samples = [group.to_numpy() for _, group in df_read.groupby(level='Tissue')]
result = f_oneway(*samples, axis=0)
df_anova = pd.DataFrame({'fstat': result.statistic, 'pval': result.pvalue}, index=df_read.columns)
df_anova.to_csv(f'{prefix}/anova.tsv', sep='\t')

# Bar chart of "significant" comparisons
alpha = 1E-10
ys = [(df_anova['pval'].dropna() < alpha).sum(),
      (df_anova['pval'].dropna() >= alpha).sum(),
      df_anova['pval'].isna().sum()]
labels = [rf'$p$ < {alpha:.0E}',
          rf'$p$ ≥ {alpha:.0E}',
          '$p$ = NaN']
xs = list(range(len(ys)))
ys, labels = list(zip(*sorted(zip(ys, labels), key=lambda x: -x[0])))  # A little bit of magic to sort them together and unzip

fig, ax = plt.subplots()
ax.bar(xs, ys, width=0.5)
ax.set_xticks(xs, labels)
ax.set_ylabel('Number of genes')
ax.set_title('$p$-values of ANOVA relative to significance level')
fig.savefig(f'{prefix}/bar|gene_number-pval.png')
plt.close()

# Histogram of F statistics
fig, ax = plt.subplots()
ax.hist(df_anova['fstat'], bins=100)
ax.set_xlabel('$F$ statistic')
ax.set_ylabel('Number of genes')
fig.savefig(f'{prefix}/hist|gene_number-Fstat.png')
plt.close()

prefix = 'out/pooled/'
if not os.path.exists(prefix):
    os.makedirs(prefix)

# Calculate t-tests (pooled)
index = df_read.index.get_level_values('Tissue') == 'liver'
sample_a = df_read[index]
sample_b = df_read[~index]
result = ttest_ind(sample_a, sample_b, axis=0, equal_var=False, alternative='greater')
df_ttest = pd.DataFrame({'tstat': result.statistic, 'pval': result.pvalue}, index=df_read.columns)
df_ttest['BH_result'] = benjamini_hochberg(df_ttest['pval'], q_BH)
df_ttest.to_csv(f'{prefix}/ttest.tsv', sep='\t')

# Pie chart of "significant" comparisons
alpha = 1E-10
ys = [(df_ttest['pval'].dropna() < alpha).sum(),
      (df_ttest['pval'].dropna() >= alpha).sum(),
      df_ttest['pval'].isna().sum()]
labels = [rf'$p$ < {alpha:.0E}',
          rf'$p$ ≥ {alpha:.0E}',
          '$p$ = NaN']
xs = list(range(len(ys)))
ys, labels = list(zip(*sorted(zip(ys, labels), key=lambda x: -x[0])))  # A little bit of magic to sort them together and unzip

fig, ax = plt.subplots()
ax.bar(xs, ys, width=0.5)
ax.set_xticks(xs, labels)
ax.set_ylabel('Number of genes')
ax.set_title('$p$-values of $t$-test relative to significance level')
fig.savefig(f'{prefix}/bar|pval.png')
plt.close()

# Histogram of p-values
fig, ax = plt.subplots()
ax.hist(df_ttest['pval'], bins=100)
ax.set_xlabel('$p$-value')
ax.set_ylabel('Number of genes')
fig.savefig(f'{prefix}/hist|gene_number-pvalue.png')
ax.set_yscale('log')
fig.savefig(f'{prefix}/hist|gene_number-pvalue|log.png')
plt.close()

# Histogram of t statistics
fig, ax = plt.subplots()
ax.hist(df_ttest['tstat'], bins=100)
ax.set_xlabel('$t$ statistic')
ax.set_ylabel('Number of genes')
fig.savefig(f'{prefix}/hist|gene_number-tstat.png')
plt.close()

# Line plot of ranked p-values
ys = df_ttest['pval'].sort_values().dropna()
xs = np.arange(1, len(ys) + 1)

fig, ax = plt.subplots()
ax.plot(xs, ys)
ax.set_xlabel('$p$-value rank')
ax.set_ylabel('$p$-value')
fig.savefig(f'{prefix}/line|pvalue-rank.png')
plt.close()

# Line plot of ranked p-values with BH cutoff
data = df_ttest[['pval', 'BH_result']].dropna().sort_values(by='pval', ignore_index=True)
index = len(data) - 1 - np.argmax(np.flip(data['BH_result'].to_numpy()))
index = int(1.05 * index)  # Expand x-axis
ys_p = data.loc[:index, 'pval']
xs_p = np.arange(1, len(ys_p) + 1)
xs_BH = np.array([1, index])
ys_BH = q_BH / len(data) * xs_BH

fig, ax = plt.subplots()
ax.plot(xs_p, ys_p, label='$p$-value')
ax.plot(xs_BH, ys_BH, label='BH critical value')
ax.set_xlabel('$p$-value rank')
ax.set_ylabel('$p$-value')
ax.legend()
fig.savefig(f'{prefix}/line|pvalue-rank|BH-rank.png')
plt.close()
