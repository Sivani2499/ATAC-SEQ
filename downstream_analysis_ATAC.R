library(tidyverse)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

# Load the preprocessed Seurat object
x.seurat <- readRDS("/Users/shivaniravindran/Desktop/single_cell_ATACSeq/x.seurat.rds")
DimPlot(x.seurat, label=TRUE)+ NoLegend()
dim(x.seurat)
#[1] 86062  7805
# 86062 peaks and 7805 cells

#Let's calculate the gene activities
gene.activities <- GeneActivity(x.seurat)
dim(gene.activities)

#add the gene activities to the seurat object as a new assay
x.seurat[['RNA']] <- CreateAssayObject(counts = gene.activities)

#normalize the gene activity matrix
x.seurat <- NormalizeData(
  object = x.seurat,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(x.seurat$nCount_RNA)
)

str(x.seurat)

#To help interpret our ATAC-seq clusters, we can visualize the activities of canonical marker genes for major cell types in PBMCs.
#Since the ATAC seq data represent measurements from sparse chromatin cell data,and because they assume a general correspondence 
#between gene body/promoter accessibility and gene expression which may not always be the case, the activities can be much noisier
#than the expression measurements from scRNA-seq data.

#Visualize the gene activities of canonical marker genes

DefaultAssay(x.seurat) <- 'RNA'
FeaturePlot(
  object = x.seurat,
  features = c('MS4A1', 'CD79A', 'CD3D', 'CD8A', 'LYZ', 'CD14', 'FCGR3A', 'NCAM1', 'CLEC4C', 'CLEC9A'),
  pt.size = 0.1,
  ncol = 5
)

#Integration with scRNA-seq data
#We can identify shared correlation patterns in the gene activity matrix and scRNA-seq dataset to identify matched biological 
#states across the two modalities. This procedure returns a classification score for each cell for each scRNA-seq-defined cluster label.
#In order to perform this integration we will be using a publicly available scRNA-seq preprocessed seurat object(based on the vignette information)

#loading the scRNA-seq seurat object
scrna <- readRDS("/Users/shivaniravindran/Downloads/pbmc_10k_v3.rds")
scrna <- UpdateSeuratObject(scrna)
str(scrna)
view(scrna@meta.data)
view(head(scrna@assays$RNA@counts))
view(scrna)

#Plotting the clusters from ATAC and RNAseq
g1 <- DimPlot(x.seurat, reduction = 'umap' )+ NoLegend() + ggtitle("scATAC-seq")
g2 <- DimPlot(scrna, reduction= 'umap', group.by = 'celltype', repel = TRUE, label = TRUE)+ NoLegend() + ggtitle("scRNA-seq")

g1 + g2

#We can see that the clusters in the scRNA-seq data are more clearly defined compared to the scATAC-seq data

#Next we will find transfer anchors between the two datasets.
transfer.anchors <- FindTransferAnchors(
  reference = scrna,
  query = x.seurat,
  reduction = 'cca'
)
#Now we will transfer the cell type labels from the scRNA-seq data to the scATAC-seq data
# This step assigns cell type labels to the ATAC-seq cells by comparing them to the labeled scRNA-seq reference.
# It uses the anchor links (found earlier) and the ATAC LSI dimensions (2–30) to find the most similar RNA cells
# and transfer their known cell type identities to the ATAC cells.
# The result is a predicted cell type label for each ATAC-seq cell.
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = scrna$celltype,
  weight.reduction = x.seurat[['lsi']],
  dims = 2:30
)
head(predicted.labels)


#Now we will add the predicted labels to the scATAC-seq object
x.seurat <- AddMetaData(object = x.seurat,metadata = predicted.labels)
#View the metadata to see the predicted labels
view(x.seurat@meta.data)

#Now we can visualize the predicted labels on the UMAP plot
p1 <- DimPlot(
  object = x.seurat,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE
) + NoLegend() + ggtitle("Predicted Cell Types from scRNA-seq")
#We can also visualize the predicted labels on the UMAP plot with the original clusters
p2 <- DimPlot(
  object = scrna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE
) + NoLegend() + ggtitle("scATAC-seq Clusters")

p1 + p2

#Let's see how many cells we have in each predicted cell type
predicted.id.counts <- table(x.seurat$predicted.id)
predicted.id.counts

#We can remove the cells that are less than 20 in number for better visualization
cells.to.keep <- names(predicted.id.counts[predicted.id.counts > 20])
cells.to.keep
#[1] "Platelet"
x.seurat <- x.seurat[,x.seurat$predicted.id %in% cells.to.keep]
#The idents look like this:GCCATAACAAACCTAC-1, so we will convert this to the predicted cell type names
Idents(x.seurat) <- x.seurat$predicted.id


#The defaultassay of the x.seurat is RNA, but we will change back to working with ATAC so seurat will focus on the ATAC peaks
DefaultAssay(x.seurat) <- 'ATAC' 
#After we’ve transferred labels and analyzed our ATAC-seq data, each cell in our x.seurat
#object now has a cell type identity (like CD4 Naive T cells, CD14+ Monocytes, etc.).
#Now we want to know which DNA regions (peaks) are more open in one cell type than in another?
# That’s what the following code does — it finds differentially accessible peaks between two groups of cells.

da.peaks <- FindMarkers(
  object = x.seurat,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14+ Monocytes",
  test.use = 'wilcox',
  min.pct = 0.1
)
head(da.peaks)
#                             p_val      avg_log2FC  pct.1. pct.2  p_val_adj
# chr2-113584614-113594917 2.640664e-240  -4.684824 0.040 0.713 2.272608e-235
# chr17-80084228-80085911  1.398974e-233   5.885278 0.399 0.010 1.203985e-228
# chr9-137263045-137268675 2.762949e-214  -4.569496 0.036 0.665 2.377850e-209
# chr7-142505185-142507751 1.379285e-209   4.841553 0.397 0.021 1.187040e-204

# The result is the table of peaks that show significant differences between the two cell types.
# Each row is a peak region, and columns mean:
# avg_log2FC — how much more open it is in one group vs the other
# 
# p_val and p_val_adj — statistical significance
# 
# pct.1 and pct.2 — fraction of cells with that peak open in each group
#For example in the first peak (chr2-113584614-113594917), the avg_log2FC is -4.684824, 
#meaning that this peak is much less accessible in CD4 Naive cells compared to CD14+ Monocytes.
#This also means that this peak is more accessible/open in CD14+ Monocytes compared to CD4 Naive cells.

#Let's take a quick peek at the distribution of the avg_log2FC values to see how many peaks are more open in each cell type
# Positive values mean the DNA region is more open in the first group that we compared (ident.1 which is CD4 naive).
# Negative values mean the DNA region is more open in the second group (ident.2 which is CD14+ monocytes).

hist(da.peaks$avg_log2FC, breaks = 100, main = "Distribution of avg_log2FC", xlab = "avg_log2FC")
#After looking at the histogram, I am setting the threshold to be >5 and <-5 for the peaks that are more open in CD4 naive and CD14 monocytes respectively
open_cd4naive <- rownames(da.peaks[da.peaks$avg_log2FC > 5, ])
open_cd14mono <- rownames(da.peaks[da.peaks$avg_log2FC < -5, ])
#Finding the closest genes to these differentially accessible peaks
closest_genes_cd4naive <- ClosestFeature(x.seurat, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(x.seurat, regions = open_cd14mono)

head(closest_genes_cd4naive)
head(closest_genes_cd14mono)

#let's visualize the top differentially accessible peaks between the two cell types using a coverage plot
top.da.peaks <- rownames(da.peaks)[1:5]
# We can also visualize an answer to the question: Do any of these open DNA regions fall near or inside the CD4 gene itself?
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd4naive), LookupGeneCoords(x.seurat, "CD4"))
CoveragePlot(
  object = x.seurat,
  region = "CD4",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)

#We can also visualize the accessibility of the LYZ gene in CD14+ Monocytes
regions_highlight <- subsetByOverlaps(StringToGRanges(open_cd14mono), LookupGeneCoords(x.seurat, "LYZ"))

CoveragePlot(
  object = x.seurat,
  region = "LYZ",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 5000
)

#we can also visualize the top 5 differentially accessible peaks between the two cell types
CoveragePlot(
  object = x.seurat,
  region = top.da.peaks,
  annotation = TRUE,
  peaks = TRUE,
  links = TRUE
)



#We can also visualize the accessibility of these peaks across all cell types using a violin plot
v1 <- VlnPlot(
  object = x.seurat,
  features = rownames(da.peaks)[1],
  pt.size = 0.1,
  group.by = 'predicted.id'
)

f1 <- FeaturePlot(
  object = x.seurat,
  features = rownames(da.peaks)[1],
  #group.by = 'predicted.id',
  pt.size = 0.1
  )
v1 + f1

#We can also create an interactive version of these plots using the CoverageBrowser() function. 

plot.ist <- CoverageBrowser(
  object = x.seurat,
  region = "CD4"
)
# This will open an interactive browser for visualizing coverage plots
