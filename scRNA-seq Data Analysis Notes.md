# Seurat Notes
## Tips
1. **Metadata** contains the general information including:   *a. Gene count b. RNA count c. Clustering ID.* Other sample information can also be stored in the metadata. Metadat can be edited by using simple assignment.
2. **Assays** contains different expression assays. Should choose the approprite assay for different purposes.



## Raw data to count matrix

**For the droplet-based method**: Each read has an unique molecular identifier (UMI) that represents the originality of the molecule (cDNA) and a cellular barcode that indicates the originality of the cell. 

**Read the folder into Seurat Object:**

```
library(Seurat)

# Read the count matrix into the Seurat object
name_file <- Read10X(data.dir = 'name')
name_object <- CreateSeuratObject(counts = name_file, min.features = 100, 
                                    min.cells = 3, project = 'name')   


# Read multiple files
for (file in c("name1", "name2")){
        seurat_data <- Read10X(data.dir = paste0("data/", file))
        seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                         min.features = 100, 
                                         project = file)
        assign(file, seurat_obj)
}

# Check each cells' RNA count 
head(name@meta.data)
```

**Merge Multiple Files to Do the Quality Control(Optional)**:
```
merged_seurat <- merge(x = name1, 
                       y = name2, 
                       add.cell.id = c("name1", "name2"))
                       
#split_seurat <- SplitObject(merged_seurat, split.by = "sample")
```



**Access the Metadata and Edit the Metadata**
```
metadata <- name_object@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)


# Add a sample column/edit column(optional)
metadata$sample <- 'sample1'
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
                
name_object@meta.data <- metadata
```

## Quality Control
**Add columns of mitochondria ratio**
```
## Add number of genes per UMI for each cell to metadata
name_object$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

## Compute percent mito ratio
name_object$mitoRatio <- PercentageFeatureSet(object = name_object, pattern = "^MT-")

name_object$mitoRatio <- name_object@meta.data$mitoRatio / 100
name_object$mitoRatio <- PercentageFeatureSet(object = name_object, pattern = "^MT-")
```

**Visulization (Optional)**
```
# General Plot
VlnPlot(pbmc, features =
c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Visualize the number of cell counts per sample
metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
metadata %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)
```

**Subset (Cell Level/Gene Level)**
```
# nUMI > 500 nGene > 250 log10GenesPerUMI > 0.8 mitoRatio < 0.2
# Filter out low quality cells using selected thresholds - these will change with experiment

#Strigent 
filtered_seurat <- subset(x = name_object, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (mitoRatio < 0.20))

#Loose
filtered_seurat <- subset(x = name_object, subset = (nFeature_RNA > 200) & (mitoRatio < 25))

# Can also subset the genes that expressed in more than 10 cells (Optional)
## Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
## Output a logical matrix specifying for each gene on whether or not there are more than zero counts per cell
nonzero <- counts > 0
## Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
## Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
## Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)


# Save the Object
save(filtered_seurat, file="data/seurat_filtered.RData")

```

## Normalization/Scale/Clustering
**Code**
```
### Individually
process_name <- NormalizeData(filtered_seurat)
process_name <- FindVariableFeatures(process_name, selection.method = 'vst',
                                nfeatures = 2000,
                                verbose = FALSE)

process_name <- ScaleData(process_name, verbose = FALSE)
process_name <- RunPCA(process_name, verbose = FALSE)

process_name <- FindNeighbors(process_name, k.param = 15, dims = 1:20)
process_name <- FindClusters(process_name, resolution = 0.5)

# Dimension Reduction
process_name <- RunTSNE(process_name, dims = 1:20)
process_name <- RunUMAP(process_name, dims = 1:20)

```
**Analysis/Visulization**
```
##### #####
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(process_name), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(process_name)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

##### To choose how many PCA dimensions to use #####
# PCA elements
print(process_name[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(process_name, dims = 1:2, reduction = "pca")


# Reduction Graph of PCA reduction
DimPlot(process_name, reduction = "pca")


# Heatmap of the PCA elements in the data
DimHeatmap(process_name, dims = 1:15, cells = 500, balanced = TRUE)

# Elbow
ElbowPlot(process_name)

##### Dimensional Plot #####
DimPlot(process_name, reduction = 'tsne', pt.size = 2, label = T, label.size = 8)

##### Find Markers #####
# 1. Individual Markers
cluster5.markers <- FindMarkers(process_name, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)

# 2. All markers
process_name.markers <- FindAllMarkers(process_name, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
process_name.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

# 3. Visulization
VlnPlot(process_name, features = c("MS4A1", "CD79A"))
FeaturePlot(process_name, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))
```

**Rename Clusters**
```
process_name <- RenameIdents(object = process_name, 
                                  "0" = "a"
                                  "1" = "b"
                                  "2" = "c")
                            
                            
## Or ##

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(process_name)
process_name <- RenameIdents(process_name, new.cluster.ids)
```
## Integration
```
# Integrate multiple seurat object

split_seurat <- list(name1, name2, name3)

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"))
 }
 
 

                                   
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
 
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
                                   
Anchors <- FindIntegrationAnchors(object.list = split_seurat)
                                            
                                            
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
```

## Reference/Helpful Links
1. [Introduction to Single-cell RNA-seq Workshop ](https://hbctraining.github.io/scRNA-seq/)
2. [Seurat Document](https://satijalab.org/seurat/)
