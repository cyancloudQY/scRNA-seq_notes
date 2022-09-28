library(Seurat)
library(tidyverse)
library(monocle3)
library(SeuratWrappers)


#### Convert seurat object to cds object
seurat_object <- readRDS('./SRR7722941_clustered_object.rds')

data <- GetAssayData(seurat_object, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat_object@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 50)     #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
plot_pc_variance_explained(cds)    #像seurat一样展示pc数

# umap
cds <- reduce_dimension(cds, preprocess_method = "PCA") #preprocess_method默认是PCA
# tSNE
cds <- reduce_dimension(cds, reduction_method="tSNE")
# 
cds <- cluster_cells(cds) 


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat_object, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed   # Make sure the UMAP is the same across Monocle and Seurat object


#Check the graph
cds <- learn_graph(cds) 
plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") 


plot_cells(cds,
           color_cells_by = "SingleR_label",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           group_label_size=4,
           cell_size=1) 





#### Plot gene expression
plot_cells(cds, color_cells_by = 'cluster', show_trajectory_graph = TRUE,  genes = 'FTL',graph_label_size=5)




### Labeled by cell type (Black Solid Circle is the Branch Nodes, Gray Circle is the Outcome)
plot_cells(cds,
           color_cells_by = "SingleR_label",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=5)


### Partition graph (whether the tree is separated)
plot_cells(cds, color_cells_by = 'partition')


########## Choose Nodes
cds <- order_cells(cds)

#### Pseudo Time
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)



######## Differential Analysis Based on Psudotime ############

cds_all_test <- graph_test(cds, neighbor_graph="principal_graph", cores=4)

deg_ids <- row.names(subset(cds_all_test, q_value < 0.05))
plot_cells(cds, color_cells_by = 'cluster', show_trajectory_graph = TRUE,  genes = 'HES4',graph_label_size=5 )


####### Find Modules of Variable Genes #######
####### Gene modules analysis (gene sets as a whole)
gene_module_df <- find_gene_modules(cds[deg_ids,], resolution=1e-2)
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(1,28)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)



########### genes
HES4_cds <- cds[rowData(cds)$gene_short_name == 'HES4',]

plot_genes_in_pseudotime(cds = HES4_cds,
                         color_cells_by="SingleR_label", min_expr = 0.05)



########## Add pseudotime to metadata
seurat_object <- AddMetaData(object = seurat_object, 
                             metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
                             col.name = 'General')
View(seurat_object@meta.data)
FeaturePlot(seurat_object, 'General') & scale_color_viridis_c()




