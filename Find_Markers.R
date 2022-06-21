library(Seurat)
library(tidyverse)
library(ggplot2)

clustered_SRR7722939 <- readRDS('./SRR7722939/data/clustered_SRR7722939.rds')

# Little Tip
# Identify Number of cells per cluster, might be useful when comparing with other sample
n_cells <- FetchData(clustered_SRR7722939, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample)  %>%
  tidyr::spread(ident, n)


# Find all markers
SRR7722939.markers <- FindAllMarkers(clustered_SRR7722939, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(SRR7722939.markers, file = './SRR7722939/results/SRR7722939_all_markers.csv')


###### Now, we need to annotate the clusters with cell type ######
## Database that contains cell markers: http://xteam.xbio.top/CellMarker/# 
## SingleRhttps://www.bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html 



#################### Visualize Markers (Manually)###########################
########## T cells ############
VlnPlot(clustered_SRR7722939, features = c("CD3D", "CD3E",'CD8A'),
        pt.size = 0.2, ncol = 3)
ggsave('./SRR7722939/results/Cell_assignment/T_cell_Vln.pdf')
FeaturePlot(clustered_SRR7722939, features = c("CD3D", "CD3E",'CD8A'),label = TRUE, pt.size = 1, label.size = 5,reduction = 'tsne')
ggsave('./SRR7722939/results/Cell_assignment/T_cell_feature.pdf')

########## B cells ############

VlnPlot(clustered_SRR7722939, features = c("CD19","CD79A","CD79B",'MS4A1'),
        pt.size = 0.2, ncol = 2)
ggsave('./SRR7722939/results/Cell_assignment/B_cell_Vln.pdf')
FeaturePlot(clustered_SRR7722939, features = c("CD19","CD79A","CD79B",'MS4A1'),label = TRUE, pt.size = 1, label.size = 5,reduction = 'tsne')
ggsave('./SRR7722939/results/Cell_assignment/B_cell_feature.pdf')

######### Monocytes & Macrophages ###########
VlnPlot(clustered_SRR7722939, features = c("CST3","LYZ","CD68"),
        pt.size = 0.2, ncol = 3)
ggsave('./SRR7722939/results/Cell_assignment/Monocytes&Macrophages_cell_Vln.pdf')
FeaturePlot(clustered_SRR7722939, features = c("CST3","LYZ","CD68"),label = TRUE, pt.size = 1, label.size = 5,reduction = 'tsne')
ggsave('./SRR7722939/results/Cell_assignment/Monocytes&Macrophages_cell_feature.pdf')


######## Fibroblasts #################

VlnPlot(clustered_SRR7722939, features = c('FGF7','MME'),
        pt.size = 0.2, ncol = 2)
ggsave('./SRR7722939/results/Cell_assignment/Fibroblasts_cell_Vln.pdf')
FeaturePlot(clustered_SRR7722939, features = c('FGF7','MME'),label = TRUE, pt.size = 1, label.size = 5,reduction = 'tsne')
ggsave('./SRR7722939/results/Cell_assignment/Fibroblasts_cell_feature.pdf')

####### Endothelial #################
VlnPlot(clustered_SRR7722939, features = c('PECAM1','VWF'),
        pt.size = 0.2, ncol = 2)
ggsave('./SRR7722939/results/Cell_assignment/Endothelial_cell_Vln.pdf')
FeaturePlot(clustered_SRR7722939, features =c('PECAM1','VWF') ,label = TRUE, pt.size = 1, label.size = 5,reduction = 'tsne')
ggsave('./SRR7722939/results/Cell_assignment/Endothelial_cell_feature.pdf')


####### NK  #########################
VlnPlot(clustered_SRR7722939, features = c('FGFBP2','FCG3RA','CX3CR1','GNLY', 'NKG7'),
        pt.size = 0.2, ncol = 2)
ggsave('./SRR7722939/results/Cell_assignment/NK_cell_Vln.pdf')
FeaturePlot(clustered_SRR7722939, features =c('FGFBP2','FCG3RA','CX3CR1','GNLY', 'NKG7') ,label = TRUE, pt.size = 1, label.size = 5,reduction = 'tsne')
ggsave('./SRR7722939/results/Cell_assignment/NK_cell_feature.pdf')

####### Tumor/Epithelial (Nah We can't.FUCK) ############
VlnPlot(clustered_SRR7722939, features = c('KRT18','KRT20'),
        pt.size = 0.2, ncol = 2)
#ggsave('./SRR7722939/results/Cell_assignment/NK_cell_Vln.pdf')
FeaturePlot(clustered_SRR7722939, features =c('KRT18','KRT20') ,label = TRUE, pt.size = 1, label.size = 5,reduction = 'tsne')
#ggsave('./SRR7722939/results/Cell_assignment/NK_cell_feature.pdf')


# It's OKAY. The original data didn't even include the tumor cell, it was peripheral blood sampel. FUCK!

#################### Let's Try SingleR to label the cells automatically ####################
library(SingleR)
library(celldex)
## load the default database
ref <- HumanPrimaryCellAtlasData()
## SingleR run
sr.run <- SingleR(test = clustered_SRR7722939@assays$RNA@data,
                  ref = ref, labels = ref$label.fine,
                  clusters = clustered_SRR7722939@active.ident,
                  fine.tune = TRUE)
## No idea what the result is but let's continue
# This is the quality control graph, need to learn later
print(plotScoreHeatmap(sr.run))

## $pruned.labels stores the cell type. Need to rename the Seurat
new.cluster.id <- sr.run$pruned.labels
names(new.cluster.id) <- levels(clustered_SRR7722939)
clustered_SRR7722939 <- RenameIdents(clustered_SRR7722939, new.cluster.id)
# We'd better add the Idents to the meta.data (for future usage)
clustered_SRR7722939@meta.data$SingleR_label <- Idents(clustered_SRR7722939)

# OK, Revisualize with the cell label
DimPlot(clustered_SRR7722939, reduction = 'umap', pt.size = 1, label = TRUE, label.size = 4)
ggsave('./SRR7722939/results/Cell_assignment/umap.cell.pdf')

DimPlot(clustered_SRR7722939, reduction = 'tsne', pt.size = 1.3, label = TRUE, label.size = 4)
ggsave('./SRR7722939/results/Cell_assignment/tsne.cell.pdf')

