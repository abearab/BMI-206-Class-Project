library(tidyverse)
library(SCENT)
library('Seurat')
library(zellkonverter)
library(SeuratDisk)
library(Matrix)
library(SingleCellExperiment)
library(GenomicRanges)
#-------------------------------------------------------------------------------
genesXcells_ALL <- readH5AD("GSE243917_genesXcells_ALL_raw.h5ad")
peaksXcells <- readH5AD("GSE243917_peaksXcells_ALL_raw.h5ad")
ATAC_ALL_CT <- readH5AD("GSE243917_ATAC_ALL_CT_processed.h5ad")
T_RNA_processed <- readH5AD("GSE243917_T_RNA_processed.h5ad")
fibroblast_processed <- readH5AD("GSE243917_fibroblast_RNA_processed.h5ad")
endothelial_processed <- readH5AD("GSE243917_endothelial_RNA_processed.h5ad")
myeloid_processed <- readH5AD("GSE243917_myeloid_RNA_processed.h5ad")
B_processed <- readH5AD("GSE243917_B_RNA_processed.h5ad")
#-------------------------------------------------------------------------------
# get RNA-seq info
rna <- assay(genesXcells_ALL, "X")
rownames(rna) <- rownames(genesXcells_ALL)  # Extract gene name
colnames(rna) <- colnames(genesXcells_ALL)
atac <- assay(peaksXcells, "X")
rownames(atac) <- rownames(peaksXcells)  # Extract peak name
colnames(atac) <- colnames(peaksXcells)  # Cell types
atac_counts <- assay(ATAC_ALL_CT, "X")

all(colnames(genesXcells_ALL) == colnames(peaksXcells))
all(colnames(peaksXcells) == colnames(ATAC_ALL_CT))
#-------------------------------------------------------------------------------
# Extract relevant columns
SCENT_obj <- CreatePeakToGeneList(SCENT_obj, genebed = "/Users/zoechen/Documents/Study/UCSF/Fall_2024/BMI_206_Biostats/Projects/GeneBody_500kb_margin.bed",
                                  nbatch = 1000,tmpfile="./temporary_atac_peak.bed",
                                  intersectedfile="./temporary_atac_peak_intersected.bed.gz")
str(SCENT_obj, max.level = 2)
# Remove rows without gene annotations
gene_peak_pairs <- gene_peak_pairs[!is.na(gene_peak_pairs$gene), ]
gene_peak_pairs <- gene_peak_pairs[gene_peak_pairs$gene %in% rownames(rna), ]
#-------------------------------------------------------------------------------
#Prepare Metadata
Tmeta <- as.data.frame(colData(T_RNA_processed))
Fmeta <- as.data.frame(colData(fibroblast_processed))
Emeta <- as.data.frame(colData(endothelial_processed))
Mmeta <- as.data.frame(colData(myeloid_processed))
Bmeta <- as.data.frame(colData(B_processed))

metadata <- rbind(Tmeta, Fmeta, Emeta, Mmeta, Bmeta)

colnames(metadata)[colnames(metadata) == "donor"] <- "sample"
colnames(metadata)[colnames(metadata) == "total_counts"] <- "nUMI"
colnames(metadata)[colnames(metadata) == "percent_mito"] <- "percent.mito"
colnames(metadata)[colnames(metadata) == "donor_num"] <- "batch"
colnames(metadata)[colnames(metadata) == "ct"] <- "celltype"
metadata$'log(nUMI)' = log(metadata$nUMI+1)
metadata$cell = rownames(metadata)

split_gene_peak_pairs <- split(gene_peak_pairs, rep(1:50, length.out = nrow(gene_peak_pairs)))

SCENT_obj <- CreateSCENTObj(rna = rna, atac = atac, meta.data = metadata, 
                            celltypes = "celltype")

SCENT_obj@peak.info.list <- split_gene_peak_pairs
#-------------------------------------------------------------------------------
# Save SCENT_obj as an RDS file
saveRDS(SCENT_obj, file = 'SCENT_500KB_nocovariates.rds')

#-------------------------------------------------------------------------------
#RUN SCENT LOCALLY
SCENT_obj <- SCENT_algorithm(object = SCENT_obj, celltype = "Tcell", ncores = 6, regr = 'poisson', bin = TRUE)

SCENT_obj@SCENT.result


