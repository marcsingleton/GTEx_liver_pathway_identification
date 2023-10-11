# GTEx Liver Pathway Analysis

This is a small analysis I did as part of a technical assessment for computational biologist role. The goal was to identify biological pathways or processes specific to the liver using the GTEx bulk RNA-seq data using only protein-coding genes and the following
tissues:

- heart
- kidney
- liver
- lung
- muscle
- pancreas
- spleen
- stomach
- pituitary gland
- thyroid

The overall approach is to first identify differentially expressed genes in the liver and then to organize these into
candidate pathways by constructing networks of co-expressed genes. Finally, term enrichment analysis of these networks
can generate hypotheses of their biological function. In more detail:

1. Data cleaning and exploratory analysis
2. Identification of differentially expressed genes in the liver
3. Construction of the co-expression network of differentially expressed genes
4. Identification of gene modules in network
5. Term enrichment analysis of genes in modules

## Project Organization

The project is organized into the following hierarchy:

```
GTEx_liver_pathway/
  ├── analysis/
  ├── bin/
  ├── data/
  │     └── GTEx/
  │           ├── bulk_count/
  │           ├── bulk_TPM/
  │           └── gene_model/
  └── src/
```