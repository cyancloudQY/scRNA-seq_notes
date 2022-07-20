library(Seurat)
library(Seurat)
library(tidyverse)
library(monocle3)
library(garnett)
library(org.Hs.eg.db)

seurat_object <- readRDS('./patch.RDS')

data <- GetAssayData(seurat_object, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat_object@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 30)
marker_check <- check_markers(cds, "./marker.csv",
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")
plot_markers(marker_check)

classifier <- train_cell_classifier(cds = cds,
                                    marker_file = "./marker.csv",
                                    db=org.Hs.eg.db,
                                    cds_gene_id_type = "SYMBOL",
                                    num_unknown = 50,
                                    marker_file_gene_id_type = "SYMBOL")

pData(cds)$garnett_cluster <- pData(cds)$seurat_clusters

cds <- classify_cells(cds, classifier, db = org.Hs.eg.db, cluster_extend = TRUE, cds_gene_id_type = "SYMBOL")
cds.meta <- subset(pData(cds), select = c("cell_type", "cluster_ext_type")) %>% as.data.frame()
View(as.data.frame(pData(cds)))


final_meta <- cbind(cell_metadata,cds.meta)
seurat_object@meta.data <- final_meta








