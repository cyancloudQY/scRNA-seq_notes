library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)
library(SingleR)
library(celldex)

## pbmc
pbmc.data <- Read10X(data.dir = "./hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2



all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")


VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
FeaturePlot(pbmc, features = c("IL7R", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
pbmc$website <- Idents(pbmc)

View(pbmc@meta.data)

Idents(pbmc) <- pbmc$seurat_clusters
######## SingleR #######
ref <- HumanPrimaryCellAtlasData()
sr.run <- SingleR(test = pbmc@assays$RNA@data,
                            ref = ref, labels = ref$label.fine,
                            clusters = pbmc@active.ident,
                            fine.tune = TRUE)

## $pruned.labels stores the cell type. Need to rename the Seurat
new.cluster.id <- sr.run$pruned.labels
names(new.cluster.id) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.id)
# We'd better add the Idents to the meta.data (for future usage)
pbmc@meta.data$SingleR_label <- Idents(pbmc)
DimPlot(pbmc, group.by = 'SingleR_label',pt.size = 1, label = TRUE)
ggsave('./Results/pbmc_singlR.pdf', width = 15, height = 10)

DimPlot(pbmc, group.by = 'website', label = TRUE, pt.size = 1)
ggsave('./Results/pbmc_website.pdf', width = 15, height = 10)

DimPlot(pbmc, group.by = 'seurat_clusters', label = TRUE, pt.size = 1)
ggsave('./Results/pbmc_raw.pdf', width = 15, height = 10)

saveRDS(pbmc, './pbmc.RDS')


########## paper1 ########  


## Clean this data matrix 
normalized_matrix <- read.table('./pancreas_refseq_rpkms_counts_3514sc.txt')
raw_count <- normalized_matrix[,c(3517:7030)]

# gene_names
all_gene_names <- normalized_matrix$V1
raw_count$gene <- all_gene_names

raw_count <-raw_count[!duplicated(raw_count$gene),]
rownames(raw_count) <- raw_count$gene
raw_count <- select(raw_count, -gene)
## Assign sample name
sample_name <- read.table('sample_name.txt')
sample_name <- t(sample_name)
colnames(raw_count) <- sample_name
seurat_object <- CreateSeuratObject(counts = raw_count)


seurat_object@meta.data$sample <- seurat_object@meta.data$orig.ident
seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA)/log10(seurat_object$nCount_RNA)
seurat_object$ERCC <- PercentageFeatureSet(object = seurat_object, pattern = "^ERCC")
seurat_object$MT <- PercentageFeatureSet(object = seurat_object, pattern = "^MT")
seurat_object$MT <- seurat_object@meta.data$MT / 100
seurat_object$Method <- 'smart'

ncells_before <- nrow(seurat_object@meta.data)

seurat_object <- subset(seurat_object, subset = ((nFeature_RNA > 2500) & (ERCC < 25) &(nCount_RNA < 200000)))


ncells_after <- nrow(seurat_object@meta.data)





seurat_object <- NormalizeData(seurat_object)
seurat_object <- ScaleData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- RunPCA(seurat_object)

ElbowPlot(seurat_object, ndims = 25)
seurat_object <- FindNeighbors(seurat_object, dims = 1: 50)
seurat_object <- FindClusters(seurat_object, resolution = 0.8)

seurat_object <- RunUMAP(seurat_object, dims = 1:50)
seurat_object <- RunTSNE(seurat_object, dims = 1:50)


DimPlot(seurat_object, reduction = 'umap', pt.size = 1, label = TRUE, label.size = 6)

saveRDS(seurat_object, './diabetes.rds')
seurat_object <- readRDS('./diabetes.rds')

# 4 acinar
VlnPlot(seurat_object, 'PRSS1')
seurat_object <- RenameIdents(seurat_object, '4' = 'acinar')

# 1 DUCTAL
VlnPlot(seurat_object, 'KRT19')
seurat_object <- RenameIdents(seurat_object, '1' = 'ductal')

# 10? MHC
VlnPlot(seurat_object, 'CD86')

# 10? mast
VlnPlot(seurat_object, 'TPSAB1')

seurat_object <- RenameIdents(seurat_object, '10' = 'MHC+MAST')

# 9 PSCs
VlnPlot(seurat_object, 'COL1A2')
seurat_object <- RenameIdents(seurat_object, '9' = 'PSCs')

# 11 endothelial 
VlnPlot(seurat_object, 'PLVAP')
seurat_object <- RenameIdents(seurat_object, '11' = 'endothelial')


# 
VlnPlot(seurat_object, 'GCG')
seurat_object <- RenameIdents(seurat_object,
                              '0' = 'alpha',
                              '2' = 'alpha',
                              '8' = 'alpha')
 
VlnPlot(seurat_object, 'INS')
seurat_object <- RenameIdents(seurat_object,
                              '5' = 'beta',
                              '7' = 'beta')


VlnPlot(seurat_object, 'PPY')
seurat_object <- RenameIdents(seurat_object,
                              '3' = 'gamma')



VlnPlot(seurat_object, 'SST')
seurat_object <- RenameIdents(seurat_object, 
                              '6' = 'delta')


VlnPlot(seurat_object, 'GHRL')

seurat_object$manual_annotate <- Idents(seurat_object)
View(seurat_object@meta.data)


saveRDS(seurat_object, './diabetes.rds')


ref <- HumanPrimaryCellAtlasData()
sr.run <- SingleR(test = seurat_object@assays$RNA@data,
                  ref = ref, labels = ref$label.fine,
                  clusters = seurat_object@active.ident,
                  fine.tune = TRUE)


new.cluster.id <- sr.run$pruned.labels
names(new.cluster.id) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.id)
# We'd better add the Idents to the meta.data (for future usage)
seurat_object@meta.data$SingleR_label <- Idents(seurat_object)

DimPlot(seurat_object, label = TRUE, pt.size = 1, group.by = 'SingleR_label')
ggsave('./Results/SingleR_diabete.pdf', width = 15, height = 10)




DimPlot(seurat_object, label = TRUE, pt.size = 1, group.by = 'manual_annotate')
ggsave('./Results/Manual_diabete.pdf', width = 12, height = 10)

DimPlot(seurat_object, label = TRUE, pt.size = 1, group.by = 'seurat_clusters')
ggsave('./Results/raw_diabete.pdf', width = 12, height = 10)

saveRDS(seurat_object, './diabetes.rds')
















