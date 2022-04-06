# you may need to install these packages; all should be available using: install.packages("package_name")
library(Seurat)
library(harmony) 
library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)
library(ggpubr)
library(RColorBrewer)

# Download AML samples: "AML419A-D0", "AML707B-D0", "AML916", and "AML921A-D0"
AML419A.D0_dem_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587950/suppl/GSM3587950_AML419A-D0.dem.txt.gz", show_col_types = FALSE)
dim(AML419A.D0_dem_tbl)
AML419A.D0_anno_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587951/suppl/GSM3587951_AML419A-D0.anno.txt.gz", show_col_types = FALSE)

AML707B.D0_dem_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587969/suppl/GSM3587969_AML707B-D0.dem.txt.gz", show_col_types = FALSE)
dim(AML707B.D0_dem_tbl)
AML707B.D0_anno_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587970/suppl/GSM3587970_AML707B-D0.anno.txt.gz", show_col_types = FALSE)

AML916.D0_dem_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587988/suppl/GSM3587988_AML916-D0.dem.txt.gz", show_col_types = FALSE)
dim(AML916.D0_dem_tbl)
AML916.D0_anno_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587989/suppl/GSM3587989_AML916-D0.anno.txt.gz", show_col_types = FALSE)

AML921A.D0_dem_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587990/suppl/GSM3587990_AML921A-D0.dem.txt.gz", show_col_types = FALSE)
dim(AML921A.D0_dem_tbl)
AML921A.D0_anno_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587991/suppl/GSM3587991_AML921A-D0.anno.txt.gz", show_col_types = FALSE)

# Download Healthy samples: "BM3", "BM4", "BM5-34p", and "BM3-34p38n"
BM3_dem_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587998/suppl/GSM3587998_BM3.dem.txt.gz", show_col_types = FALSE)
dim(BM3_dem_tbl)
BM3_anno_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3587nnn/GSM3587999/suppl/GSM3587999_BM3.anno.txt.gz", show_col_types = FALSE)

BM4_dem_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3588nnn/GSM3588000/suppl/GSM3588000_BM4.dem.txt.gz", show_col_types = FALSE)
dim(BM4_dem_tbl)
BM4_anno_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3588nnn/GSM3588001/suppl/GSM3588001_BM4.anno.txt.gz", show_col_types = FALSE)

BM5.34p_dem_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3588nnn/GSM3588002/suppl/GSM3588002_BM5-34p.dem.txt.gz", show_col_types = FALSE)
dim(BM5.34p_dem_tbl)
BM5.34p_anno_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3588nnn/GSM3588002/suppl/GSM3588002_BM5-34p.anno.txt.gz", show_col_types = FALSE)

BM5.34p38n_dem_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3588nnn/GSM3588003/suppl/GSM3588003_BM5-34p38n.dem.txt.gz", show_col_types = FALSE)
dim(BM5.34p38n_dem_tbl)
BM5.34p38n_anno_tbl = read_tsv("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3588nnn/GSM3588003/suppl/GSM3588003_BM5-34p38n.anno.txt.gz", show_col_types = FALSE)

# process data frames
AML419A.D0_dem_tbl <- AML419A.D0_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()
AML707B.D0_dem_tbl <- AML707B.D0_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()
AML916.D0_dem_tbl <- AML916.D0_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()
AML921A.D0_dem_tbl <- AML921A.D0_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()

BM3_dem_tbl <- BM3_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()
BM4_dem_tbl <- BM4_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()
BM5.34p_dem_tbl <- BM5.34p_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()
BM5.34p38n_dem_tbl <- BM5.34p38n_dem_tbl %>% column_to_rownames("Gene") %>% as.matrix()

AML419A.D0_anno_tbl <- AML419A.D0_anno_tbl %>% column_to_rownames("Cell") 
AML707B.D0_anno_tbl <- AML707B.D0_anno_tbl %>% column_to_rownames("Cell")
AML916.D0_anno_tbl <- AML916.D0_anno_tbl %>% column_to_rownames("Cell") 
AML921A.D0_anno_tbl <- AML921A.D0_anno_tbl %>% column_to_rownames("Cell")

BM3_anno_tbl <- BM3_anno_tbl %>% column_to_rownames("Cell")
BM4_anno_tbl <- BM4_anno_tbl %>% column_to_rownames("Cell")
BM5.34p_anno_tbl <- BM5.34p_anno_tbl %>% column_to_rownames("Cell")
BM5.34p38n_anno_tbl <- BM5.34p38n_anno_tbl %>% column_to_rownames("Cell")

# create seurat objects for each sample
aml1 = CreateSeuratObject(counts = AML419A.D0_dem_tbl, min.cells = 1, min.features = 1)
aml2 = CreateSeuratObject(counts = AML707B.D0_dem_tbl, min.cells = 1, min.features = 1)
aml3 = CreateSeuratObject(counts = AML916.D0_dem_tbl, min.cells = 1, min.features = 1)
aml4 = CreateSeuratObject(counts = AML921A.D0_dem_tbl, min.cells = 1, min.features = 1)

bm1 = CreateSeuratObject(counts = BM3_dem_tbl, min.cells = 1, min.features = 1)
bm2 = CreateSeuratObject(counts = BM4_dem_tbl, min.cells = 1, min.features = 1)
bm3 = CreateSeuratObject(counts = BM5.34p_dem_tbl, min.cells = 1, min.features = 1)
bm4 = CreateSeuratObject(counts = BM5.34p38n_dem_tbl, min.cells = 1, min.features = 1)

# add Sample name and group in metadata of each seurat object
aml1$Sample <- "AML419A-D0"
aml2$Sample <- "AML707B-D0"
aml3$Sample <- "AML916-D0"
aml4$Sample <- "AML921A-D0"

bm1$Sample <- "BM3"
bm2$Sample <- "BM4"
bm3$Sample <- "BM5-34p"
bm4$Sample <- "BM5-34p38n"

aml1$Group <- "AML"
aml2$Group <- "AML"
aml3$Group <- "AML"
aml4$Group <- "AML"

bm1$Group <- "Healthy"
bm2$Group <- "Healthy"
bm3$Group <- "Healthy"
bm4$Group <- "Healthy"

# combine Seurat objects into one object
combo <- merge(aml1, c(aml2, aml3, aml4, bm1, bm2, bm3, bm4))

# create metadata dataframe with all sample metadata
dim(AML419A.D0_anno_tbl) # has an extra column for NanoporeTranscripts (we'll need to remove this to merge the metadata data frames)
dim(AML707B.D0_anno_tbl) # has an extra column for NanoporeTranscripts (we'll need to remove this to merge the metadata data frames)
dim(AML916.D0_anno_tbl)
dim(AML921A.D0_anno_tbl)
dim(BM3_anno_tbl)
dim(BM4_anno_tbl)
dim(BM5.34p_anno_tbl)
dim(BM5.34p38n_anno_tbl)

AML419A.D0_anno_tbl <- AML419A.D0_anno_tbl[,1:27] # remove the extra column from this data frame
AML707B.D0_anno_tbl <- AML707B.D0_anno_tbl[,1:27] # remove the extra column from this data frame

# combine all the meta data sheets into one data frame
metadata <- rbind(AML419A.D0_anno_tbl, AML707B.D0_anno_tbl, AML916.D0_anno_tbl, AML921A.D0_anno_tbl,
      BM3_anno_tbl, BM4_anno_tbl, BM5.34p_anno_tbl, BM5.34p38n_anno_tbl)

# order the metadata dataframe so that it matches the barcode order in the seurat object
metadata <- metadata[match(colnames(combo), rownames(metadata)),]

# add all columns from metadata to seurat object
for(x in colnames(metadata)){
  combo[[x]] <- metadata[,x]
}

## add percentage of reads mapping to mitochondrial genome in metadata column
combo[["percent.mito"]] = PercentageFeatureSet(combo, pattern = "^MT-")

# plot qc metrics:
p1 <- VlnPlot(combo, features = "nFeature_RNA", pt.size = 0) + 
  labs(x = "Sample", y = "Number of genes detected per cell", title = "") + 
  scale_fill_brewer(palette = "Paired") + NoLegend() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2 <- VlnPlot(combo, features = "nCount_RNA", pt.size = 0) + 
  labs(x = "Sample", y = "Number of UMIs per cell", title = "") + 
  scale_fill_brewer(palette = "Paired") + NoLegend() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p3 <- VlnPlot(combo, features = "percent.mito", pt.size = 0) + 
  labs(x = "Sample", y = "Percentage of reads mapping\nto the mitochondrial genome", title = "") + 
  scale_fill_brewer(palette = "Paired") + NoLegend() + 
  scale_y_continuous(labels = function(x) paste0(x, "%"))

tiff("qc_plots.tiff", units = "in", height = 12, width = 6, res = 300) # note if you're not familiar: the tiff() function will save the plot as a .tiff file to your working directory
plot_grid(p1, p2, p3, nrow = 3, rel_heights = c(1, 1, 1.4))
dev.off()
 
# To combine samples and remove technical differences between them, we perform integration with Harmony
# Igor didn't show us how to integrate between samples; here I'm using an integration tool called Harmony. Here's a vignette that describes it more: https://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html
# there's also an integration vignette on the Seurat website which does this slightly differently, but it takes longer

# first we find variable features; we want to find variable features in each sample individually, and then choose the top 3000 variable genes that are found among the most samples
# the SelectIntegrationFeatures function does this on a list of seurat objects; so we first split our data into a list by Sample
object.list <- SplitObject(combo, split.by = "Sample")
features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = 3000) # find top 3000 variable genes (can use 2 or 3k)

# Log-normalize data
combo <- NormalizeData(combo) 

# scale and regress out technical variables (helps with creating a well-mixed UMAP, but not always necessary); note that we scale our top variable features we just identified using SelectIntegrationFeatures
combo <- ScaleData(combo, features = features, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito")) 

# run PCA, again using our variable features
combo <- RunPCA(combo, features = features, npcs = 50, verbose = TRUE)

# Integrate samples using Harmony
set.seed(12) # make sure to set a random seed so integration results are 100% reproducible (otherwise the results are stochastic)
combo <- combo %>% RunHarmony("Sample", plot_convergence = TRUE)

# not super important, but this PCA elbow plot helps us choose the number of PCs to select for integration, and shows how these PCs are corrected for batch differences, if you compare the two plots (unintegrated vs integrated)
tiff("integration_corrected_pca_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(combo, ndims = 50, reduction = "harmony")
dev.off()

tiff("uncorrected_pca_elbowplot.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(combo, ndims = 50, reduction = "pca") 
dev.off()

# we'll chooose 20 PCs based on the elbow plot for dimensionality reduction with UMAP
combo <- RunUMAP(combo, reduction = "harmony", dims = 1:20)

# color the UMAP by Group (AML vs Healthy), Sample, and Malignency, as well as the CellType the authors assigned to each cell
DimPlot(combo, group.by = "Group", reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(combo, group.by = "Sample", reduction = "umap", cols = RColorBrewer::brewer.pal(8, "Paired")) + theme(aspect.ratio = 1)
DimPlot(combo, group.by = "PredictionRefined", reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(combo, group.by = "CellType", reduction = "umap") + theme(aspect.ratio = 1)

# highlight each of these cells on the UMAP (not necessary, just if you're interested)
Idents(combo) <- combo$CellType
for(x in setdiff(unique(combo$CellType), NA)){
  cells <- WhichCells(combo, idents = x)
  p <- DimPlot(combo, cells.highlight = cells) + ggtitle(x) + theme(aspect.ratio = 1)
  plot(p)
}

# to compare our integrated results with the results we'd get if we didn't integrate the data, we also run UMAP algorithm on our non-batch corrected PCA results
combo <- RunUMAP(combo, reduction = "pca", dims = 1:20, reduction.name = "umap_unintegrated") # save un-integrated UMAP results under reductions as "umap_unintegrated"
DimPlot(combo, group.by = "Sample", reduction = "umap_unintegrated", cols = RColorBrewer::brewer.pal(8, "Paired")) + theme(aspect.ratio = 1)
DimPlot(combo, group.by = "Group", reduction = "umap_unintegrated") + theme(aspect.ratio = 1)
DimPlot(combo, group.by = "CellType", reduction = "umap_unintegrated") + theme(aspect.ratio = 1)

## directly compare umap results from unintegrated versus integrated UMAPs

# first by sample; these results will show how integration helps us find similar cell types across different samples (integrated UMAP is much more mixed across samples)
leg <- get_legend(DimPlot(combo, group.by = "Sample", reduction = "umap_unintegrated", cols = RColorBrewer::brewer.pal(8, "Paired")) + theme(aspect.ratio = 1) + labs(x = "UMAP1", y = "UMAP2", title = "Unintegrated")) # get legend 
p1 <- DimPlot(combo, group.by = "Sample", reduction = "umap_unintegrated", cols = RColorBrewer::brewer.pal(8, "Paired")) + theme(aspect.ratio = 1) + labs(x = "UMAP1", y = "UMAP2", title = "Unintegrated") + NoLegend() # get unintegrated plot
p2 <- DimPlot(combo, group.by = "Sample", reduction = "umap", cols = RColorBrewer::brewer.pal(8, "Paired")) + theme(aspect.ratio = 1) + labs(x = "UMAP1", y = "UMAP2", title = "Integrated") + NoLegend() # get integrated plot
tiff("umap_integrated_vs_unintegrated_bysample.tiff", units = "in", width = 15, height = 6, res = 300) # combine them and save to your working directory
plot_grid(p1, p2, leg, rel_widths = c(1, 1, 0.5), nrow = 1)
dev.off()

# make pretty color palette
tol21rainbow=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

# same comparison plot but with author's identified cell types
leg <- get_legend(DimPlot(combo, group.by = "CellType", reduction = "umap_unintegrated", cols = tol21rainbow) + theme(aspect.ratio = 1) + labs(x = "UMAP1", y = "UMAP2", title = "Unintegrated"))
p1 <- DimPlot(combo, group.by = "CellType", reduction = "umap_unintegrated", cols = tol21rainbow) + theme(aspect.ratio = 1) + labs(x = "UMAP1", y = "UMAP2", title = "Unintegrated") + NoLegend()
p2 <- DimPlot(combo, group.by = "CellType", reduction = "umap", cols = tol21rainbow) + theme(aspect.ratio = 1) + labs(x = "UMAP1", y = "UMAP2", title = "Integrated") + NoLegend()
tiff("umap_integrated_vs_unintegrated_bycelltype.tiff", units = "in", width = 15, height = 6, res = 300)
plot_grid(p1, p2, leg, rel_widths = c(1, 1, 0.5), nrow = 1)
dev.off()

# Create nearest neighbor graph (for clustering); note we use same number of PCs as when we run UMAP
combo <- FindNeighbors(combo, reduction = "harmony", dims = 1:20)

# run K-means clustering algorithm (I included several resolutions, which control how many clusters to make)
combo <- FindClusters(combo, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))

# I looked at all the resolutions but 0.5 seemed to give the most similar results to the cell types in the paper; feel free to check the others; they're saved in the Seurat object metadata as RNA_snn_res.<resolution>
Idents(combo) <- combo$RNA_snn_res.0.5
DimPlot(combo, reduction = "umap", cols = tol21rainbow) # plot UMAP with clustering results

# the authors assigned cells predicted as malignant to be "-like"; 
#for our purposes, let's combine them, so we can see our broad cell type clusters
DimPlot(combo, reduction = "umap", cols = tol21rainbow, group.by = "CellType")

Idents(combo) <- combo$CellType # temporarily make CellType the main ident
test = RenameIdents(combo,"cDC-like"="cDC", # assign new idents (I made a new seurat object called "test" for this purpose)
                        "GMP-like"="GMP",
                        "HSC-like"="HSC",
                        "Mono-like"="Mono",
                        "Prog-like"="Prog",
                        "ProMono-like"="ProMono")

DimPlot(test, reduction = "umap", cols = tol21rainbow)
DimPlot(test, reduction = "umap", cells.highlight = colnames(test[,is.na(Idents(test))])) # weirdly, one cell is unlabeled

combo$paper_idents <- test@active.ident # assign these new identifiers to a metadata column in the original seurat object called "paper_idents"

# from the metadata spreadsheets we downloaded for each sample, the authors assigned each cell a score for a number of different cell types
# these scores are derived from expression of gene sets canonically associated with each cell type
# we can test to see what our cell types our clusters might belong to using these gene set scores

# visualize all the gene sets (and my preliminary notes)
VlnPlot(combo, features = "CyclingScore", cols = tol21rainbow) # 0 and 1
VlnPlot(combo, features = "Score_HSC", cols = tol21rainbow) # 0 and 1
VlnPlot(combo, features = "Score_Prog", cols = tol21rainbow) # 0, 1, (slightly) 4/5, 7, 8, 11
VlnPlot(combo, features = "Score_GMP", cols = tol21rainbow) # 8
VlnPlot(combo, features = "Score_ProMono", cols = tol21rainbow) # 3
VlnPlot(combo, features = "Score_Mono", cols = tol21rainbow) # 2
VlnPlot(combo, features = "Score_cDC", cols = tol21rainbow) # 4
VlnPlot(combo, features = "Score_pDC", cols = tol21rainbow) # 4 [need to subcluster to identify cDC vs pDC]
VlnPlot(combo, features = "Score_earlyEry", cols = tol21rainbow) # 5
VlnPlot(combo, features = "Score_lateEry", cols = tol21rainbow) # 5 (may need to further split to get early v late)
VlnPlot(combo, features = "Score_ProB", cols = tol21rainbow) # 11
VlnPlot(combo, features = "Score_B", cols = tol21rainbow) # 10
VlnPlot(combo, features = "Score_Plasma", cols = tol21rainbow) # 12
VlnPlot(combo, features = "Score_T", cols = tol21rainbow) # 6
VlnPlot(combo, features = "Score_CTL", cols = tol21rainbow) # 6 & 9 [probably 9]
VlnPlot(combo, features = "Score_NK", cols = tol21rainbow) # 9 [9 appears to be a mix of CTL and NK cells]

# likely cluster assignments
# Cluster 0 = HSC (because high in HSC and Prog genes but not cycling)
# Cluster 1 = Prog (because high in HSC and Prog genes and cycling)
# Cluster 2 = Monocytes
# Cluster 3 = Promonocytes (because high in mono and promono genes and cycling)
# Cluster 4 = cDC
# Cluster 5 = Erythrocyte (may need to further cluster to split into early and late)
# Cluster 6 = T
# Cluster 7 = Prog
# Cluster 8 = GMP (high in Prog genes, and GMP genes)
# Cluster 9 = CTL & NK (need to further cluster to split)
# Cluster 10 = B cells
# Cluster 11 = Pro-B (proB genes and cycling)
# Cluster 12 = Plasma

# missing early/late erythrocyte divide
# missing NK cells (merged with T-cells right now
# the missing clusters were likely not resolved at our clustering resolution because these are smaller cell type lineage differences;
# to pull out these missing cell types, I subclustered the cell types likely to contain these missing clusters

# subclustering requires us to subset the Seurat object so we have a new seurat object containing just the cells in the cluster we wish to subcluster
# we'll subcluster cluster 9 (to get 2 clusters); and subcluster cluster 5 (to get 2 clusters)

# let's start with cluster 5, which appears to contain erythrocytes, but doesn't have a clear early/late divide in our original clustering
cluster5 <- subset(combo, idents = 5)

# annoyingly, if we want to pull out differences in an individual cluster, we need to start over identifying new variable features, scaling, running PCA, and integrating
# these steps are the same as we did earlier, but for just our individual cluster seurat object

# first we find sample specific variable features
cluster.5.object.list <- SplitObject(cluster5, split.by = "Sample")
cluster5.features <- SelectIntegrationFeatures(object.list = cluster.5.object.list, nfeatures = 2000)

# scale and regress out technical variables
cluster5 <- ScaleData(cluster5, features = cluster5.features, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito")) # regressing out percent.mt 

# run PCA
cluster5 <- RunPCA(cluster5, features = cluster5.features, npcs = 50, verbose = TRUE)

# Integrate samples using Harmony
set.seed(333) # make sure to set a random seed so integration results are 100% reproducible
cluster5 <- cluster5 %>% RunHarmony("Sample", plot_convergence = TRUE)

# make cluster 5 subclustering elbow plot; you'll notice, most variability is found in only the first few PCs; we'll use 6 PCs for downstream analysis
tiff("integration_corrected_pca_elbowplot_cluster5.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(cluster5, ndims = 50, reduction = "harmony")
dev.off()

tiff("uncorrected_pca_elbowplot_cluster5.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(cluster5, ndims = 50, reduction = "pca") 
dev.off()

# run UMAP using first 6 PCs
cluster5 <- RunUMAP(cluster5, reduction = "harmony", dims = 1:6)

# color UMAP by features
DimPlot(cluster5, group.by = "Group", reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(cluster5, group.by = "Sample", reduction = "umap", cols = RColorBrewer::brewer.pal(8, "Paired")) + theme(aspect.ratio = 1)
DimPlot(cluster5, group.by = "PredictionRefined", reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(cluster5, group.by = "CellType", reduction = "umap", cols = tol21rainbow) + theme(aspect.ratio = 1) # you'll see our UMAP separates early/late Erys pretty well

# run clustering 
cluster5 <- FindNeighbors(cluster5, reduction = "harmony", dims = 1:6)
cluster5 <- FindClusters(cluster5, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))

# I chose 0.3 as the resolution, since I want to get only 2 clusters, early and late Erys
Idents(cluster5) <- cluster5$RNA_snn_res.0.3
DimPlot(cluster5, reduction = "umap") # plot UMAP results

# make violin plots for classifying these subclusters
VlnPlot(cluster5, features = "Score_earlyEry", cols = tol21rainbow) # Cluster 1 = early erythro
VlnPlot(cluster5, features = "Score_lateEry", cols = tol21rainbow) # Cluster 0 is very clearly enriched in late ery genes

# rename these clusters with their likely idents
cluster5 = RenameIdents(cluster5, '0' = "lateEry",
                        '1'="earlyEry")

### next do the same thing with Cluster 9; which looks to be a mix of CTLs and NKs
cluster9 <- subset(combo, idents = 9)

# first we find sample specific variable features
cluster.9.object.list <- SplitObject(cluster9, split.by = "Sample")
cluster9.features <- SelectIntegrationFeatures(object.list = cluster.9.object.list, nfeatures = 2000)

# scale and regress out technical variables
cluster9 <- ScaleData(cluster9, features = cluster9.features, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito")) # regressing out percent.mt 

# run PCA
cluster9 <- RunPCA(cluster9, features = cluster9.features, npcs = 50, verbose = TRUE)

# Integrate samples using Harmony
set.seed(936) # make sure to set a random seed so integration results are 100% reproducible
cluster9 <- cluster9 %>% RunHarmony("Sample", plot_convergence = TRUE)

# cluster 9 elbow plots; we'll use first 8 PCs for UMAP
tiff("integration_corrected_pca_elbowplot_cluster9.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(cluster9, ndims = 50, reduction = "harmony")
dev.off()

tiff("uncorrected_pca_elbowplot_cluster9.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(cluster9, ndims = 50, reduction = "pca") 
dev.off()

# run UMAP
cluster9 <- RunUMAP(cluster9, reduction = "harmony", dims = 1:8)

# color UMAP by cell types; we'll see CTLs and NKs are pretty well separated by our UMAP
DimPlot(cluster9, group.by = "Group", reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(cluster9, group.by = "Sample", reduction = "umap", cols = RColorBrewer::brewer.pal(8, "Paired")) + theme(aspect.ratio = 1)
DimPlot(cluster9, group.by = "PredictionRefined", reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(cluster9, group.by = "CellType", reduction = "umap", cols = tol21rainbow) + theme(aspect.ratio = 1)

# run clustering
cluster9 <- FindNeighbors(cluster9, reduction = "harmony", dims = 1:8)
cluster9 <- FindClusters(cluster9, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))

Idents(cluster9) <- cluster9$RNA_snn_res.0.3 # split into 5 clusters; splitting less leaves some likely NKs out
DimPlot(cluster9, reduction = "umap")

# create violin plots of gene set scores, to help us assign our clusters to cell types
VlnPlot(cluster9, features = "Score_T", cols = tol21rainbow) 
VlnPlot(cluster9, features = "Score_CTL", cols = tol21rainbow) # Cluster 0 & 3 are clearly CTL
VlnPlot(cluster9, features = "Score_NK", cols = tol21rainbow) # Cluster 1 & 2 are clearly NK

# cluster 4 is not high in either of these scores; for now let's call it Unknown

# we'll assign Clusters 1 & 2 as NK, 0 & 3 as CTL, and label Cluster 4 as "Unknown" for now
cluster9 = RenameIdents(cluster9, '0' = "CTL",
                    '1'="NK",
                    '2'="NK",
                    '3'="CTL",
                    '4'="Unknown")

# visualize our final classifications; they look pretty good
VlnPlot(cluster9, features = "Score_T", cols = tol21rainbow) 
VlnPlot(cluster9, features = "Score_CTL", cols = tol21rainbow) 
VlnPlot(cluster9, features = "Score_NK", cols = tol21rainbow)
VlnPlot(cluster9, features = "Score_Prog", cols = tol21rainbow)

#### also subcluster cluster 4:
### next do the same thing with Cluster 9; which looks to be a mix of CTLs and NKs
cluster4 <- subset(combo, idents = 4)

# first we find sample specific variable features
cluster.4.object.list <- SplitObject(cluster4, split.by = "Sample")
cluster4.features <- SelectIntegrationFeatures(object.list = cluster.4.object.list, nfeatures = 2000)

# scale and regress out technical variables
cluster4 <- ScaleData(cluster4, features = cluster4.features, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mito")) # regressing out percent.mt 

# run PCA
cluster4 <- RunPCA(cluster4, features = cluster4.features, npcs = 50, verbose = TRUE)

# Integrate samples using Harmony
set.seed(936) # make sure to set a random seed so integration results are 100% reproducible
cluster4 <- cluster4 %>% RunHarmony("Sample", plot_convergence = TRUE)

# cluster 9 elbow plots; we'll use first 10 PCs for UMAP
tiff("integration_corrected_pca_elbowplot_cluster4.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(cluster4, ndims = 50, reduction = "harmony")
dev.off()

tiff("uncorrected_pca_elbowplot_cluster4.tiff", units = "in", width = 5, height = 5, res = 300)
ElbowPlot(cluster4, ndims = 50, reduction = "pca") 
dev.off()

# run UMAP
cluster4 <- RunUMAP(cluster4, reduction = "harmony", dims = 1:10)

# color UMAP by cell types; we'll see CTLs and NKs are pretty well separated by our UMAP
DimPlot(cluster4, group.by = "Group", reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(cluster4, group.by = "Sample", reduction = "umap", cols = RColorBrewer::brewer.pal(8, "Paired")) + theme(aspect.ratio = 1)
DimPlot(cluster4, group.by = "PredictionRefined", reduction = "umap") + theme(aspect.ratio = 1)
DimPlot(cluster4, group.by = "CellType", reduction = "umap", cols = tol21rainbow) + theme(aspect.ratio = 1)

# run clustering
cluster4 <- FindNeighbors(cluster4, reduction = "harmony", dims = 1:10)
cluster4 <- FindClusters(cluster4, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))

Idents(cluster4) <- cluster4$RNA_snn_res.0.1 # split into 3 clusters # 1 and 2 appear to be cDC, 0 is pDC
DimPlot(cluster4, reduction = "umap")

# create violin plots of gene set scores, to help us assign our clusters to cell types
VlnPlot(cluster4, features = "Score_cDC", cols = tol21rainbow) # Cluster 1 & 2
VlnPlot(cluster4, features = "Score_pDC", cols = tol21rainbow) # Cluster 0 (though mean differences are not huge)

# cluster 4 is not high in either of these scores; for now let's call it Unknown

# we'll assign Clusters 1 & 2 as NK, 0 & 3 as CTL, and label Cluster 4 as "Unknown" for now
cluster4 = RenameIdents(cluster4, '0' = "pDC",
                        '1'="cDC",
                        '2'="cDC")

# visualize our final classifications; they look pretty good
VlnPlot(cluster4, features = "Score_cDC", cols = tol21rainbow)
VlnPlot(cluster4, features = "Score_pDC", cols = tol21rainbow)


### now we'll need to add these identifiers back to the main Seurat object

# first get lists of cells in each final cell type cluster
NK.cells <- WhichCells(cluster9, idents = "NK")
CTL.cells <- WhichCells(cluster9, idents = "CTL")
Unknown.cells <- WhichCells(cluster9, idents = "Unknown")
earlyEry.cells <- WhichCells(cluster5, idents = "earlyEry")
lateEry.cells <- WhichCells(cluster5, idents = "lateEry")
pDC.cells <- WhichCells(cluster4, idents = "pDC")
cDC.cells <- WhichCells(cluster4, idents = "cDC")

# make a dataframe of the main seurat objects clustering results
cluster.ident.df <- data.frame("barcode" = colnames(combo), 
                               "cluster" = as.character(combo@active.ident))

# we'll then rename the cells which are found in our lists of new identifiers
cluster.ident.df <- cluster.ident.df %>% mutate(modified = ifelse(barcode %in% NK.cells, "NK", 
                                                    ifelse(barcode %in% CTL.cells, "CTL",
                                                           ifelse(barcode %in% Unknown.cells, "Unknown",
                                                                  ifelse(barcode %in% earlyEry.cells, "earlyEry",
                                                                         ifelse(barcode %in% lateEry.cells, "lateEry",
                                                                                ifelse(barcode %in% pDC.cells, "pDC",
                                                                                       ifelse(barcode %in% cDC.cells, "cDC",
                                                                                cluster))))))))
# finally, assign new cluster identifiers to main seurat object under metadata field "modified"
combo$modified <- cluster.ident.df$modified

# rename remaining clusters
# Cluster 0 = HSC (because high in HSC and Prog genes but not cycling)
# Cluster 1 = Prog (because high in HSC and Prog genes and cycling)
# Cluster 2 = Monocytes
# Cluster 3 = Promonocytes (because high in mono and promono genes and cycling)
# Cluster 4 = cDC
# Cluster 5 = Erythrocyte (may need to further cluster to split into early and late)
# Cluster 6 = T
# Cluster 7 = Prog
# Cluster 8 = GMP (high in Prog genes, and GMP genes)
# Cluster 9 = CTL & NK (need to further cluster to split)
# Cluster 10 = B cells
# Cluster 11 = Pro-B (proB genes and cycling)
# Cluster 12 = Plasma
Idents(combo) <- combo$modified

combo = RenameIdents(combo, '0' = "HSC",
                        '1'="Prog",
                        '2'="Mono",
                        '3'="ProMono",
                        '6'="T",
                        '7'="Prog",
                        '8'="GMP",
                        '10'="B",
                        '11'="ProB",
                        '12'="Plasma")

# save these final identifiers to "final_ident" metadata column
combo$final_ident <- combo@active.ident

# visualize results on UMAP plot and compare to original paper_idents; our clusters look more cohesive, a good sign
DimPlot(combo, reduction = "umap", cols = tol21rainbow)
DimPlot(combo, reduction = "umap", cols = tol21rainbow, group.by = "paper_idents")

# highlight the unknown cluster to examine where it is on UMAP
cells <- WhichCells(combo, idents = "Unknown")
DimPlot(combo, reduction = "umap", cells.highlight = cells)

# to figure out what this unknown cluster is, we perform DE testing between Unknown cells and all other cells
unknown.markers <- FindMarkers(combo, ident.1 = "Unknown")
unknown.markers <- unknown.markers %>% arrange(desc(avg_log2FC)) # several granzyme genes; these are very high in CTLs; this is a CTL cluster

# rename Unknown cells as CTL
combo = RenameIdents(combo, 'Unknown'="CTL")
combo$final_ident <- combo@active.ident # save final idents

## setting levels so the order stays the same on our plots
combo$paper_idents <- factor(combo$paper_idents, levels = c("B", "cDC", "CTL", "earlyEry", "GMP", "HSC", "lateEry", "Mono", "NK", "pDC", "Plasma", "ProB", "Prog", "ProMono", "T"))
DimPlot(combo, group.by = "paper_idents", cols = tol21rainbow)

combo$final_ident <- factor(combo$final_ident, levels = c("B", "cDC", "CTL", "earlyEry", "GMP", "HSC", "lateEry", "Mono", "NK", "pDC", "Plasma", "ProB", "Prog", "ProMono", "T"))
combo@active.ident <- factor(combo@active.ident, levels = c("B", "cDC", "CTL", "earlyEry", "GMP", "HSC", "lateEry", "Mono", "NK", "pDC", "Plasma", "ProB", "Prog", "ProMono", "T"))

# plot unintegrated results for author assigned idents, & UMAP results
leg = get_legend(DimPlot(combo, group.by = "paper_idents", reduction = "umap_unintegrated", cols = tol21rainbow))
p1 = DimPlot(combo, group.by = "paper_idents", reduction = "umap_unintegrated", cols = tol21rainbow) + ggtitle("Unintegrated") + labs(x = "UMAP_1", y = "UMAP_2") + theme(aspect.ratio = 1) + NoLegend()
p2 = DimPlot(combo, group.by = "paper_idents", cols = tol21rainbow) + ggtitle("Integrated") + theme(aspect.ratio = 1) + NoLegend()

tiff("final_paperidents_by_integration.tiff", units = "in", width = 14, height = 6, res = 300)
plot_grid(p1, p2, leg, nrow = 1, rel_widths = c(1,1,0.2))
dev.off()

# plot our final cluster assignments on UMAPs, both unintegrated and integrated (For comparison)
leg = get_legend(DimPlot(combo, reduction = "umap_unintegrated", cols = tol21rainbow))
p1 = DimPlot(combo, reduction = "umap_unintegrated", cols = tol21rainbow) + ggtitle("Unintegrated") + labs(x = "UMAP_1", y = "UMAP_2") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + NoLegend()
p2 = DimPlot(combo, cols = tol21rainbow) + ggtitle("Integrated") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + NoLegend()

tiff("final_clusters_by_integration.tiff", units = "in", width = 14, height = 6, res = 300)
plot_grid(p1, p2, leg, nrow = 1, rel_widths = c(1,1,0.2))
dev.off()

# plot malignence on UMAPs, both unintegrated and integrated (For comparison)
leg = get_legend(DimPlot(combo, reduction = "umap_unintegrated", group.by = "PredictionRefined", cols = c("red3", "skyblue")))
p1 = DimPlot(combo, reduction = "umap_unintegrated", group.by = "PredictionRefined", cols = c("red3", "skyblue")) + ggtitle("Unintegrated") + labs(x = "UMAP_1", y = "UMAP_2") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + NoLegend()
p2 = DimPlot(combo, group.by = "PredictionRefined", cols = c("red3", "skyblue")) + ggtitle("Integrated") + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) + NoLegend()

tiff("malignence_by_integration.tiff", units = "in", width = 14, height = 6, res = 300)
plot_grid(p1, p2, leg, nrow = 1, rel_widths = c(1,1,0.2))
dev.off()

# plot original author assigned cell types
tiff("umap_celltypes_authors.tiff", units = "in", width = 8, height = 8, res = 300)
DimPlot(combo, group.by = "CellType", cols = tol21rainbow) + ggtitle("Cell Types (by authors)") +  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
dev.off()

# plot original author assigned cell types (grouping malignent cells with healthy)
tiff("umap_celltypes_authors_malig_norm_combined.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(combo, group.by = "paper_idents", cols = tol21rainbow) + ggtitle("Cell Types (by authors)", subtitle = "normal + malignant") +  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5),
                                                                                                                                            plot.subtitle= element_text(hjust = 0.5, face="italic"))
dev.off()

# plot our final clusters
tiff("umap_ourclusters.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(combo, group.by = "final_ident", cols = tol21rainbow) + ggtitle("Cell Types (our clusters)") +  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
dev.off()

# plot samples
tiff("umap_sample.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(combo, group.by = "Sample", cols = RColorBrewer::brewer.pal(8, "Paired")) + ggtitle("Sample") +  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
dev.off()

# plot UMAP with malignence
tiff("umap_malignence.tiff", units = "in", width = 6, height = 6, res = 300)
DimPlot(combo, group.by = "PredictionRefined", cols = c("red3", "skyblue")) + ggtitle("Malignance") +  theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
dev.off()


### let's create violin plots to validate our clustering
VlnPlot(combo, features = "CyclingScore", cols = tol21rainbow, pt.size = 0) #
VlnPlot(combo, features = "Score_HSC", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_Prog", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_GMP", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_ProMono", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_Mono", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_cDC", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_pDC", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_earlyEry", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_lateEry", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_ProB", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_B", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_Plasma", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_T", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_CTL", cols = tol21rainbow, pt.size = 0) # 
VlnPlot(combo, features = "Score_NK", cols = tol21rainbow, pt.size = 0) # 

# compare to violin plot score results for the author's original cell type identifiers (note these are very similar to ours)
VlnPlot(combo, features = "CyclingScore", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) #
VlnPlot(combo, features = "Score_HSC", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_Prog", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_GMP", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_ProMono", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_Mono", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_cDC", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_pDC", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) #
VlnPlot(combo, features = "Score_earlyEry", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_lateEry", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_ProB", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_B", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_Plasma", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_T", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_CTL", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 
VlnPlot(combo, features = "Score_NK", cols = tol21rainbow, group.by = "paper_idents", pt.size = 0) # 

# now we'll create boxplots to compare our results, since these are a little nicer than violins for summarizing stats

# create datafraem with score results based on our cluster assignments
tmp.df <- data.frame("Cluster" = combo@active.ident,
                     "Score_HSC" = combo$Score_HSC,
                     "Score_Prog" = combo$Score_Prog,
                     "Score_GMP" = combo$Score_GMP,
                     "Score_ProMono" = combo$Score_ProMono,
                     "Score_Mono" = combo$Score_Mono,
                     "Score_cDC" = combo$Score_cDC,
                     "Score_pDC" = combo$Score_pDC,
                     "Score_earlyEry" = combo$Score_earlyEry,
                     "Score_lateEry" = combo$Score_lateEry,
                     "Score_ProB" = combo$Score_ProB,
                     "Score_B" = combo$Score_B,
                     "Score_Plasma" = combo$Score_Plasma,
                     "Score_T" = combo$Score_T,
                     "Score_CTL" = combo$Score_CTL,
                     "Score_NK" = combo$Score_NK,
                     "Analysis" = rep("Ours", ncol(combo)))

# same thing but with author's original assignments
tmp.df.2 <- data.frame("Cluster" = combo$paper_idents, 
                     "Score_HSC" = combo$Score_HSC,
                     "Score_Prog" = combo$Score_Prog,
                     "Score_GMP" = combo$Score_GMP,
                     "Score_ProMono" = combo$Score_ProMono,
                     "Score_Mono" = combo$Score_Mono,
                     "Score_cDC" = combo$Score_cDC,
                     "Score_pDC" = combo$Score_pDC,
                     "Score_earlyEry" = combo$Score_earlyEry,
                     "Score_lateEry" = combo$Score_lateEry,
                     "Score_ProB" = combo$Score_ProB,
                     "Score_B" = combo$Score_B,
                     "Score_Plasma" = combo$Score_Plasma,
                     "Score_T" = combo$Score_T,
                     "Score_CTL" = combo$Score_CTL,
                     "Score_NK" = combo$Score_NK,
                     "Analysis" = rep("Theirs", ncol(combo)))

# combine them for plotting
tmp.df.3 <- rbind(tmp.df, tmp.df.2)

# filter our the Cluster authors assigned as NA originally
tmp.df.3 <- tmp.df.3 %>% filter(!is.na(Cluster)) # exclude NA cells

# for loop saving plots for each of these cluster scores; 
for(x in levels(combo@active.ident)){
  tiff(paste0(x, "_score_boxplot_analysis.tiff"), units = "in", height = 6, width = 6, res = 300)
  plot(ggboxplot(tmp.df.3, x = "Cluster", y = paste0("Score_", x),
            color = "black", fill = "Analysis", notch = TRUE, palette=c("purple3", "orange"),
            outlier.shape = NA) + labs(y = paste0(x, " Score"), title = x) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) + coord_flip()
  )
  dev.off()
}

saveRDS(combo, "combo_final_seurat_object.rds")

##### Differential expression testing (for 5 DE tests will likely be necessary)
## 1. differential expression testing across our clusters
cluster.markers <- FindAllMarkers(combo) # run DE testing
cluster.markers <- cluster.markers %>% arrange(cluster, desc(avg_log2FC)) # sort by cluster and log2 fold change 
# cluster.markers %>% filter(cluster == "B") %>% top_n(10, avg_log2FC) # example of pulling out top 10 B cell cluster markers

## 2. differential expression testing across the authors clusters
Idents(combo) <- combo$CellType
authors.celltypes.markers <- FindAllMarkers(combo)
authors.celltypes.markers <- authors.celltypes.markers %>% arrange(cluster, desc(avg_log2FC))

## 3. differential expression testing across the authors clusters (combining malignant + normal cells)[grouping -like's with cell type clusters]
Idents(combo) <- combo$paper_idents
authors.celltypes.combined.malig.norm.markers <- FindAllMarkers(combo)
authors.celltypes.combined.malig.norm.markers <- authors.celltypes.combined.malig.norm.markers %>% arrange(cluster, desc(avg_log2FC))

## 4. differential expression testing within each cluster between malignant and healthy cells
# haven't done this yet, but it'll look something like this: 
# Idents(combo) <- combo$final_ident
# for(x in unique(Idents(combo))){
#   markers <- FindMarkers(pbmc_small, ident.1 = "g1", group.by = 'groups', subset.ident = "x")
# }

## 5. differential expression testing within each cell type (from the original paper) between malignant and healthy cells
# haven't done this yet

## ****haven't gotten to this part yet: 
## Heatmaps showing top DEGs between those 5 DE tests
##
## 1. differential expression testing across our clusters

## 2. differential expression testing across the authors clusters


## 3. differential expression testing across the authors clusters (combining malignant + normal cells)


## 4. differential expression testing within each cluster between malignant and healthy cells


## 5. differential expression testing within each cell type (from the original paper) between malignant and healthy cells







