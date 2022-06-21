library(Seurat)
library(tidyverse)
library(ggplot2)


# Read in the filtered file
filtered_SRR7722939 <- readRDS('./SRR7722939/data/filtered_SRR7722939.rds')

############### Normalizing the data ###########################
## Normalized data stored in """filtered_SRR7722939[['RNA']]@data"""
## By the way, counts are stored in """filtered_SRR7722939[['RNA']]@counts"""
filtered_SRR7722939 <- NormalizeData(filtered_SRR7722939)

## Variable Feature
filtered_SRR7722939 <- FindVariableFeatures(filtered_SRR7722939, selection.method = 'vst', nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(filtered_SRR7722939), 10)

# plot variable features with and without labels
# Use """FetchData(filtered_SRR7722939, vars = c('S100A9'))""" to get the expression of the variable genes
plot1 <- VariableFeaturePlot(filtered_SRR7722939)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
ggsave('./SRR7722939/results/top10_variable.pdf')



################ Scaling the data #############################
#1. The scale data can be used to do the dimensional reduction
#2. Can not use this Scaled data to do the expression analysis or 
#visulization, because they are scaled to around o and 1. 
#3. Use normalzied data for analysis.
filtered_SRR7722939 <- ScaleData(filtered_SRR7722939, vars.to.regress = 'mitoRatio')


################### Run PCA on Data ###########################
# 1. To reduce the dimension
# 2. PCA scores served as a meta feature for the later analysis
# 3. Gotta choose how many PCs to use

filtered_SRR7722939 <- RunPCA(filtered_SRR7722939, features = VariableFeatures(object = filtered_SRR7722939))
### Some Visualization
print(filtered_SRR7722939[['pca']], dims = 1:5, nfeatures = 5)

# PCA1,PCA2 genes
VizDimLoadings(filtered_SRR7722939, dims = 1:2, reduction = "pca")
ggsave('./SRR7722939/results/pc1_pc2_genes_rank.pdf')
# Just a dimensional reduction plot
DimPlot(filtered_SRR7722939, reduction = "pca")
ggsave('./SRR7722939/results/pc1_vs_pc2_red.pdf')
# Heatmap (might be a little bit useful for correlation analysis by serving as a reference)
# Generally, same color means same trend in these genes
pdf('./SRR7722939/results/pca_heatmap.pdf')
DimHeatmap(filtered_SRR7722939, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

## Still hard to choose how many PCs to use for later analysis? 
## Don't worry. Do JackStraw Procedure and draw a elbow plot. 
filtered_SRR7722939 <- JackStraw(filtered_SRR7722939, num.replicate = 100)

filtered_SRR7722939 <- ScoreJackStraw(filtered_SRR7722939, dims = 1:20)
JackStrawPlot(filtered_SRR7722939, dims = 1:15)
ggsave('./SRR7722939/results/jack_straw.pdf')


ElbowPlot(filtered_SRR7722939)

ggsave('./SRR7722939/results/elbow_plot.pdf')

######################################################################
## So 12 looks like a good cut off based on these PCA graphs :] ######
######################################################################
###################### Cluster the Cells #############################
### 1. During the Following process, I found 12 PCs bring too much noise
### 2. So I chose 8 PCs to rerun and get an ideal UMAP graph
### 3. Always remember these parameters can be tuned for the best visualization

# graph-based algorithm
# Oh, don't forget to use 12 dimensions we determined from PCA
# Screw it, 8 seems like a better number after a test run

filtered_SRR7722939 <- FindNeighbors(filtered_SRR7722939,dims = 1:8)

# But it's really hard to choose the resolution right here
# emmmmm. Let's choose 0.8 here first. 
filtered_SRR7722939 <- FindClusters(filtered_SRR7722939, resolution = 0.8)
# Use """Idents(filtered_SRR7722939)""" to check the cell Identity, aka CLUSTER
# OK, let's try 1 and 1.2
filtered_SRR7722939 <- FindClusters(filtered_SRR7722939, resolution = 1)
filtered_SRR7722939 <- FindClusters(filtered_SRR7722939, resolution = 1.2)

## Run UMAP/tSNE for visualization purpose. Definitely remember to use the dimensions we chose!!!
filtered_SRR7722939 <- RunUMAP(filtered_SRR7722939, dims = 1:8)
filtered_SRR7722939 <- RunTSNE(filtered_SRR7722939, dims = 1:8)



## So smooth if we just use 8 dimensions instead of 12.
# We can always change the Idents of the cells for visualization purpose/cluster purpose
# Change Idents won't change the graph's view, because graph is based on UMAP. 
Idents(filtered_SRR7722939)  <- 'RNA_snn_res.1.2'
DimPlot(filtered_SRR7722939, reduction = 'umap', pt.size = 1, label = TRUE, label.size = 6)
ggsave('./SRR7722939/results/UMAP_raw.pdf')


# Don't know why, but TSNE plot looks more smoooooooth.Might Use TSNE plots in downstream analysis
DimPlot(filtered_SRR7722939, reduction = 'tsne', pt.size = 1.5, label = TRUE, label.size = 6)
ggsave('./SRR7722939/results/TSNE_raw.pdf')


saveRDS(filtered_SRR7722939, file = "./SRR7722939/data/clustered_SRR7722939.rds")



