library(Seurat)
library(ggplot2)
library(tidyverse)
library(DoubletFinder)

## 42
counts <- Read10X('./filtered_feature_bc_matrix//')
seurat_object <- CreateSeuratObject(counts = counts, min.features = 100)
seurat_object$mitoRatio <- PercentageFeatureSet(object = seurat_object, pattern = "^MT-")
seurat_object$mitoRatio <- seurat_object@meta.data$mitoRatio / 100
VlnPlot(seurat_object, features = c('nFeature_RNA', 'nCount_RNA', 'mitoRatio'), ncol = 3, group.by =  'orig.ident')

before_cn <- ncol(seurat_object)

seurat_object <- subset(seurat_object, subset = (nFeature_RNA > 200) &
                          (nCount_RNA > 500) &
                          (mitoRatio < 0.25))

after_cn <- ncol(seurat_object)


seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
seurat_object <- ScaleData(seurat_object)
seurat_object <- RunPCA(seurat_object)
ElbowPlot(seurat_object)
seurat_object <- RunUMAP(seurat_object, dims = 1:10)

seurat_object <- FindNeighbors(seurat_object, dims = 1:15)
seurat_object <- FindClusters(seurat_object, resolution = 1)
DimPlot(seurat_object)

########### Doublet Finder ##############
sweep.res.list <- paramSweep_v3(seurat_object, PCs = 1:15, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)

mpk<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
# mpk = 0.05

## Estimate Doublet Number
#### Estimate homotypic doublet
annotations <- seurat_object$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 

#### Every 1000 cells there is gonna be a 8/1000 increasing  
# Doublet Rate = number of cells * 8 * 1e-6 
doublet_rate <- ncol(seurat_object) * 8 * 1e-6 
nExp_poi <- doublet_rate * ncol(seurat_object)
nExp_poi.adj <- nExp_poi * (1 - homotypic.prop)



##### Do it 
seurat_object <- doubletFinder_v3(seurat_object, PCs = 1:15, pN = 0.25, pK = mpk, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
seurat_object <- doubletFinder_v3(seurat_object, PCs = 1:15, pN = 0.25, pK = mpk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
View(seurat_object@meta.data)



######## Scrublet
scrublet <- read.csv('./doublet_prediction.txt', header = T, sep = ',')
scrublet$prediction[scrublet$prediction == 'False']  <- 'Singlet'
scrublet$prediction[scrublet$prediction == 'True'] <- 'Doublet'


meta.data <- seurat_object@meta.data
meta.data$barcode <- rownames(meta.data)
meta.data <- left_join(meta.data, scrublet, by = c('barcode' = 'barcodes'))
rownames(meta.data) <- meta.data$barcode

seurat_object@meta.data <- meta.data


DimPlot(seurat_object, group.by = 'DF.classifications_0.25_0.005_119.482656', pt.size = 1) + labs(title = 'Doublet Finder Results')
ggsave('./Results/Doublet_Finder_UMAP.pdf', width = 10, height = 8) 

DimPlot(seurat_object, group.by = 'prediction' , pt.size = 1) + labs(title = 'Scrublet Results')
ggsave('./Results/Scrublet_UMAP.pdf', width = 10, height = 8) 




seurat_object$all_results <- 'Singlet'
seurat_object$all_results[seurat_object$DF.classifications_0.25_0.005_119.482656 == 'Doublet'] <- 'DF'
seurat_object$all_results[seurat_object$prediction == 'Doublet'] <- 'SCR'
seurat_object$all_results[which(seurat_object$prediction == 'Doublet' & seurat_object$DF.classifications_0.25_0.005_119.482656 == 'Doublet')] <- 'Both'


DimPlot(seurat_object, group.by = 'all_results', cols = c('red','gold','cyan','black'), pt.size = 1) +
  labs(title = 'Integrated Graph')
ggsave('./Results/Both_UMAP.pdf', width = 10, height = 8) 


table(seurat_object$prediction)
table(seurat_object$DF.classifications_0.25_0.005_119.482656)




###############################
df <- subset(seurat_object, subset = DF.classifications_0.25_0.005_119.482656 == 'Doublet')
scr <- subset(seurat_object, subset = prediction == 'Doublet')

df <- rownames(df@meta.data)
scr <- rownames(scr@meta.data)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    DF = df,
    Scrublet = scr
  ),
  filename = "Venn.png", imagetype = 'png',
  col = "transparent",
  fill = c("red", "blue"),
  alpha = 0.5,
  label.col = c("darkred", "white", 'red'),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.default.pos = "text",
  cat.col = c("darkred", "darkblue"),
  cat.cex = 2.5,
  cat.fontfamily = "serif",
  cat.dist = c(0.06, 0.06),
  cat.pos = 0
)

########### Some analysis  ###########
seurat_object$general_results <- 'Singlet'
seurat_object$general_results[which(seurat_object$all_results %in% c('DF', 'Both', 'SCR'))] <- 'Doublet'
View(seurat_object@meta.data)
table(seurat_object$general_results)

VlnPlot(seurat_object,'nFeature_RNA',group.by = 'general_results') + 
  geom_boxplot(width=0.2,position=position_dodge(0.9)) + 
  stat_boxplot(geom="errorbar",width=0.1,size=0.5,position=position_dodge(0.6),color="black") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave('./Results/Violin_Gene_Number.pdf', width = 6, heigh = 5)

VlnPlot(seurat_object,'nCount_RNA',group.by = 'general_results') + 
  geom_boxplot(width=0.2,position=position_dodge(0.9)) + 
  stat_boxplot(geom="errorbar",width=0.1,size=0.5,position=position_dodge(0.6),color="black") + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave('./Results/Violin_Transcript_Number.pdf', width = 6, heigh = 5)


########## 假设检验 ###########
sta_test <- select(seurat_object@meta.data, c('nFeature_RNA', 'general_results'))

# Shapiro-Wilk normality test (Normality?)
with(sta_test, shapiro.test(nFeature_RNA[general_results == "Singlet"]))
with(sta_test, shapiro.test(nFeature_RNA[general_results == "Doublet"]))

# deviation (similar?)
res.ftest <- var.test(nFeature_RNA ~ general_results, data = sta_test)
res.ftest


##### Can't use t-test. Use wilcox-test
res <- wilcox.test(nFeature_RNA ~ general_results, data = sta_test)
res # p-value < 2.2e-16

##### Count
sta_test_2 <- select(seurat_object@meta.data, c('nCount_RNA', 'general_results'))
res <- wilcox.test(nCount_RNA ~ general_results, data = sta_test_2)
res # p-value < 2.2e-16


###### SingleR ##########
library(celldex)
library(SingleR)

ref <- HumanPrimaryCellAtlasData()
sr.run <- SingleR(test = seurat_object@assays$RNA@data,
                  ref = ref, labels = ref$label.fine,
                  clusters = seurat_object@active.ident,
                  fine.tune = TRUE)

## $pruned.labels stores the cell type. Need to rename the Seurat
new.cluster.id <- sr.run$pruned.labels
names(new.cluster.id) <- levels(seurat_object)
seurat_object <- RenameIdents(seurat_object, new.cluster.id)
# We'd better add the Idents to the meta.data (for future usage)
seurat_object@meta.data$SingleR_label <- Idents(seurat_object)
DimPlot(seurat_object, group.by = 'SingleR_label',pt.size = 0.7, label = TRUE)
ggsave('./SingleR.pdf', width = 13, height = 10)







