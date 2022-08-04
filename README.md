<!-- README.md is generated from README.Rmd. Please edit that file -->

# Data Analysis on Metabolomics Data

## :book: Procedures

### Statistical Analysis

Due to terrible experience on *Statistical Analysis in Metabolomics* via
**MetaboAnalystR** R package, we try to provide a reproducible and
easy-to-use template for visualization, pre-processing, exploration, and
statistical analysis on metabolomic data by other packages and scripts.
Here, the template comprises the following procedures:

1.  Data Processing

    -   Data Checking

    -   Data Filtering

    -   Missing Value Imputation

    -   Data Normalization

2.  Cluster Analysis

    -   Hierarchical Clustering

    -   Partitional Clustering

3.  Chemometrics Analysis

    -   Principal Component Analysis (PCA)

    -   Partial Least Squares-Discriminant Analysis (PLS-DA)

    -   Sparse Partial Least Squares-Discriminant Analysis (sPLS-DA)

4.  Univariate Analysis

    -   Fold Change Analysis

    -   T Tests

    -   Wilcoxon Test

    -   Limma Test

    -   Wilcoxon Test

    -   Volcano plot

    -   Correlation Heatmaps

    -   glasso

5.  Feature selection

    -   Lasso

    -   Ridge

    -   Elasticnet

6.  Classification

    -   Random Forest

7.  Network Analysis

    -   SPRING

    -   Spearman

    -   SparCC

    -   Network comparison

### Functional Analysis

**Following two chapters would focus on the Enrichment Analysis and
Pathway Analysis of metabolomic data. Enrichment Analysis includes three
sections (i.e., ORA, SSP and QEA) and Pathway Analysis only includes ORA
and QEA.**

**The main difference between Enrichment Analysis and Pathway Analysis
are the data set that input metabolites are enriched to. In Enrichment
Analysis, input metabolites are enriched to pre-defined metabolite sets
while in Pathway Analysis, metabolites are enriched to pathways in
KEGG.**

1.  Enrichment Analysis

    -   Single Sample Profiling

    -   Over representation analysis

    -   Quantitative Enrichment Analysis

2.  Pathway Analysis

    -   Over representation analysis

    -   Quantitative Enrichment Analysis

## :writing_hand: Authors

1.  [Hua Zou](zouhua@xbiome.com)

2.  [Bangzhuo Tong](zouhua@xbiome.com)

Xbiome company

## :wrench: Change log

-   Submitted to gitlab. (2022-06-28)
-   add `test.Rmd`. (2022-07-05)
-   add `template.Rmd`. (2022-07-08)
-   add *README.Rmd*. (2022-07-12)
-   add *building bookdwon*. (2022-08-04)
