library(DESeq2)
library(rafalib) #BiocManager::install("rafalib")
library(ggplot2) #BiocManager::install("ggplot2")
library(pheatmap) #BiocManager::install("pheatmap")
library(RColorBrewer) #BiocManager::install("RColorBrewer")
library(dplyr) #BiocManager::install("dplyr")
library(plotly) #BiocManager::install("plotly")
library(org.Hs.eg.db) #BiocManager::install("org.Hs.eg.db")
library(devtools) #BiocManager::install("devtools")
library(EnsDb.Hsapiens.v86) #BiocManager::install("EnsDb.Hsapiens.v86")
library(gridExtra) #BiocManager::install("gridExtra")
library(EnhancedVolcano) #BiocManager::install("EnhancedVolcano")
library(survminer) #BiocManager::install("survminer")
library(biomaRt) #BiocManager::install("biomaRt")
library(tidyverse) #BiocManager::install("tidyverse")
library(genefilter) #BiocManager::install("genefilter")
library(tidyr) #BiocManager::install("tidyr")
library(VennDiagram) #BiocManager::install("VennDiagram")
library(ggbeeswarm)

setwd()

#INPUT data
countdata_main <- read.table(
    "raw_counts.txt", 
    header=TRUE, row.names=1
    )

countdata_mirna <- read.table(
    "miRNA_raw_counts.txt", 
    header=TRUE, row.names=1
    )
condition_main <- read.csv(
    "sample_metadata_final.csv",
    header=TRUE, row.names=1, sep=","
    )

##Making dds element
(coldata_main <- data.frame(
    row.names=colnames(countdata_main), condition_main))

##Feed this info into the dds matrix element
dds_main <- DESeqDataSetFromMatrix(
    countData = countdata_main, 
    colData = coldata_main, 
    design = ~ Cohort
    )

dds_main <- estimateSizeFactors(dds_main)

#Variant stabilizing transformation of the data
vsd_main <- varianceStabilizingTransformation(dds_main, blind=FALSE)
head(assay(vsd_main))
hist(assay(vsd_main))


dds_main <- DESeq(dds_main)


#print normalized counts for the upload
dds_main_counts<-counts(dds_main, normalized=TRUE)
write.table(dds_main_counts, file="dds_main_counts.txt")



res.LFCshrink <- lfcShrink(dds_main, contrast=c("Cohort", "IBM", "Control_Amputee"), type = "ashr")

write.table(res.LFCshrink, "res.LFCshrink.txt")

#PCA plot

mat <- assay(vsd_main)
assay(vsd_main) <- mat

#PCA plots using DESeq2 plot function
##Cohort
PCA_Cohort <- plotPCA(vsd_main, intgroup = c("Cohort"), returnData = T)
PCA_Cohort <- plotPCA(vsd_main, intgroup = c("Cohort"))
show(PCA_Cohort)
ggplot(PCA_Cohort, aes(x=PC1,y=PC2,col=Cohort,label=name)) 






#NO NEED TO RUN THE CODE BELOW


#renaming the fist column to GeneID and remove ensembl version from the identifiers
LFCmain_new <- tibble::rownames_to_column(LFCmain, "GeneID")
#convertiing the txt file to csv
write.csv(LFCmain_new, "LFCmain_new.csv", quote = F)
LFCmain_new$GeneID <- sub('(?<=\\.).*$', '', LFCmain_new$GeneID, perl=TRUE)
LFCmain_new$GeneID <- substr(LFCmain_new$GeneID, 1, nchar(LFCmain_new$GeneID)-1)
#row.names(LFCmain_new) <- LFCmain_new$GeneID
#LFCmain_new[1] <- NULL
# keep only the significant genes
results_sig = subset(LFCmain_new, padj < 0.01)
#write.csv(results_sig, "significant.csv")





