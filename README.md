# Master-Thesis-Inclusion-Body-Myositis
In this repository you can find the scripts I used for my Master Thesis in Inclusion Body Myositis

# Run the analysis, by starting from the differential gene expression, followed by the enrichment analysis, the WGCNA and finally the eQTL analysis

# Differential Gene Expression analysis
Simple differential gene expression analysis using DESeq2. This script was produced by another member within the ejprd group.

# GO enrichment analysis
Enriched GO terms of the differentially expressed genes

# HGNC symbols
Add HGNC symbols next to the ensembl identifiers, using the biomart package

# Weighted gene co-expression analysis (WGCNA)
A script to find the co-expressed genes in IBM

# eQTL analysis
To run this analysis, you will need the whole exome sequencing (WES) dataset, that is currently not available publicly. Therefore, the analysis
cannot be done without this dataset.



# Information about the files

raw_counts.txt = raw counts of the samples used for the analysis

res.LFCshrink.txt = results of the differentiall gene expression analysis, using the log-fold shrinkage argument

sample metadata.txt = covariates of the samples

normalized_id_eQTLs_input.txt = contains the normalized genes identifiers, and their chromosomes starting/ending points

#For any question, please send an email at t.zarotiadis@student.maastrichtuniversity.nl


