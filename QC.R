library(Seurat)
library(tidyverse)
library(ggplot2)
####### Read in the count matrix and create Seurat Object
counts <- Read10X(data.dir = './SRR7722939/raw_feature_bc_matrix/')
SRR7722939 <- CreateSeuratObject(counts = counts, min.features = 100)
# Maybe add a sample column?
SRR7722939$sample <- 'SRR7722939'


######### QULITY CONTROL ############
## Add some columns
# Number of genes per UMI
SRR7722939$log10GenesPerUMI <- log10(SRR7722939$nFeature_RNA)/log10(SRR7722939$nCount_RNA)

# Mitochondria genes ratio
SRR7722939$mitoRatio <- PercentageFeatureSet(object = SRR7722939, pattern = "^MT-")
SRR7722939$mitoRatio <- SRR7722939@meta.data$mitoRatio / 100

## Some Quality Control Plots as long as you want 
## (One of them can be used to show the filtering cutoff))
# Summary plot
VlnPlot(SRR7722939, features = c('nFeature_RNA', 'nCount_RNA', 'mitoRatio'), ncol = 3)
ggsave('./SRR7722939/QC_results/num of genes & mitoRatio.pdf', width = 10, height = 4)

# num of cells
SRR7722939@meta.data %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells") +
  geom_text(aes(label=..count..),stat = 'count')
ggsave('./SRR7722939/QC_results/num of cells.pdf')


## Can be used to drewn quality control figure
# Visualize the number UMIs/transcripts per cell
plot1 <- SRR7722939@meta.data %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
plot2 <- SRR7722939@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

# Visualize the distribution of mitochondrial gene expression detected per cell
plot3 <- SRR7722939@meta.data %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# Visualize nFeature RNA
plot4 <- SRR7722939@meta.data %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 200)

plot1 + plot2 + plot3 +plot4
ggsave('./SRR7722939/QC_results/countRNA&complexity&mitoRatio.pdf')


## Reeaaaaal Subset Step
# Stringent filtering (need further consideration :>)
# Cell Level
#filtered_SRR7722939 <- subset(x = SRR7722939, subset = ((nCount_RNA >= 500)&
#                                                              (nFeature_RNA >= 250) &
#                                                              (log10GenesPerUMI > 0.8) &
#                                                              (mitoRatio < 0.2)))

# Gene Level (Should we even do this step??? Maybe not if we don't care about average expression of a cell)
#counts <- GetAssayData(object = filtered_SRR7722939, slot = "counts")
#nonzero <- counts > 0
#keep_genes <- Matrix::rowSums(nonzero) >= 10
#filtered_counts <- counts[keep_genes, ]
#filtered_SRR7722939 <- CreateSeuratObject(filtered_counts, meta.data = filtered_SRR7722939@meta.data)


# Loose filtering (Then we keep more cells :|) 
filtered_SRR7722939  <- subset(x = SRR7722939, subset = (nFeature_RNA > 200) & (mitoRatio < 0.25))
VlnPlot(filtered_SRR7722939, features = c('nFeature_RNA', 'nCount_RNA', 'mitoRatio'), ncol = 3)
ggsave('./SRR7722939/QC_results/after_filtering.pdf', width = 10, height = 4)

# Cell number after

filtered_SRR7722939@meta.data %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells") +
  geom_text(aes(label=..count..),stat = 'count')
ggsave('./SRR7722939/QC_results/num of cells_after.pdf')


saveRDS(filtered_SRR7722939, file="./SRR7722939/data/filtered_SRR7722939.rds")





