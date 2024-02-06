"""Plot statistics related to identification of differentially expressed genes in the liver."""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from src.constants import q_BH

parser = argparse.ArgumentParser()
parser.add_argument('anova')
parser.add_argument('ttest')
parser.add_argument('corr')
parser.add_argument('-o', '--output_path', default='./')

if __name__ == '__main__':
      args = parser.parse_args()

      df_anova = pd.read_table(args.anova)
      df_ttest = pd.read_table(args.ttest)
      df_corr = pd.read_table(args.corr, header=[0, 1], index_col=[0, 1])

      prefix = f'{args.output_path}/unpooled/'
      if not os.path.exists(prefix):
            os.makedirs(prefix)

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

      prefix = f'{args.output_path}/pooled/'
      if not os.path.exists(prefix):
            os.makedirs(prefix)

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

      # Histogram of correlations
      array = df_corr.to_numpy().ravel()

      fig, ax = plt.subplots()
      ax.hist(array, bins=100)
      ax.set_xlabel('Correlation coefficient')
      ax.set_ylabel('Number of gene pairs')
      fig.savefig(f'{prefix}/hist|pair_number-corr.png')
      plt.close()
