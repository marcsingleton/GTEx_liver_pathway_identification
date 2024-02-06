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

import argparse
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, ttest_ind
from src.constants import tissue_regex, q_BH


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


parser = argparse.ArgumentParser()
parser.add_argument('gtf_path')
parser.add_argument('read_path')
parser.add_argument('-o', '--output_path', default='./')

if __name__ == '__main__':
    args = parser.parse_args()

    # Load gene metadata
    df_gene = pd.read_table(args.gtf_path)
    genes_protein = df_gene.loc[df_gene['gene_type'] == 'protein_coding', 'gene_id']
    genes_other = df_gene.loc[df_gene['gene_type'] != 'protein_coding', 'gene_id']

    # Load mapped reads
    dfs_read = []
    files = sorted([file for file in os.listdir(args.read_path) if file.endswith('.gct.gz')])
    for file in files:
        # Read in data and transpose to make genes columns
        df_read = pd.read_table(f'{args.read_path}/{file}', skiprows=2)
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

    prefix = f'{args.output_path}/unpooled/'
    if not os.path.exists(prefix):
        os.makedirs(prefix)

    # Calculate ANOVAs (unpooled)
    samples = [group.to_numpy() for _, group in df_read.groupby(level='Tissue')]
    result = f_oneway(*samples, axis=0)
    df_anova = pd.DataFrame({'fstat': result.statistic, 'pval': result.pvalue}, index=df_read.columns)
    df_anova.to_csv(f'{prefix}/anova.tsv', sep='\t')

    prefix = f'{args.output_path}/pooled/'
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

    # Subset to liver genes and calculate correlations
    genes_liver = df_ttest[df_ttest['BH_result']].index.get_level_values('Name')
    df_liver = df_read.loc[df_read.index.get_level_values('Tissue') == 'liver', genes_liver]
    df_corr = df_liver.corr()
    df_corr.to_csv(f'{prefix}/corr.tsv', sep='\t')
