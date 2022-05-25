library(Seurat)
library(Matrix)
data <- readRDS(snakemake@input[[1]])
DefaultAssay(object = data) <- "RNA"
write.csv(rownames(data),snakemake@output[['gene_names']])
write.csv(data@reductions$Spatial@cell.embeddings,snakemake@output[['location']])
writeMM(data@assays$RNA@counts,file=snakemake@output[['expression']])