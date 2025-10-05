# script to process single-cell ATAC-Seq data
# Vignette: https://stuartlab.org/signac/articles/pbmc_vignette
getwd()
setwd("/Users/shivaniravindran/Desktop/single_cell_ATACSeq/")

# install packages
# remotes::install_github("stuart-lab/signac", ref="develop")
BiocManager::install("EnsDb.Hsapiens.v75")
BiocManager::install("biovizBase")

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(tidyverse)
library(biovizBase)

#reading the fragment file
fragmen <- read.delim('/Users/shivaniravindran/Downloads/atac_pbmc_10k_v1_fragments.tsv.gz', header = F, nrows = 10)
head(fragmen)

# 1. Read in the .h5 file

counts <- Read10X_h5('/Users/shivaniravindran/Downloads/atac_pbmc_10k_v1_raw_peak_bc_matrix.h5')
counts[1:10,1:10]

#creating a chromatin assay object 

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "/Users/shivaniravindran/Downloads/atac_pbmc_10k_v1_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

str(chrom_assay)

#reading the metadat - this is usually the cell barcode metrics data and will be in csv format
metadata <- read.csv(file = '/Users/shivaniravindran/Downloads/atac_pbmc_10k_v1_singlecell.csv', header = T, row.names = 1)
View(metadata)

#create seurat obje
x.seurat <- CreateSeuratObject(
 counts = chrom_assay,
 meta.data = metadata,
 assay = "ATAC"
)

str(x.seurat)

#checking for the annotation slot - you will see this empty for now
x.seurat@assays$ATAC@annotation

##extract the gene annotation form th ensembl database. 
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
annotations

#when you view the annotations you will see that the seqnames are in a different formt
#Inorder to change the format of this to the UCSC, we will be adding the string 'str' for each row in the seqnames column

seqlevels(annotations) <- paste0("chr", seqlevels(annotations))

#Now we will add this annotations object to the x.seurat object that we created earlier using the function Annotation()

Annotation(x.seurat) <- annotations
x.seurat@assays$ATAC@annotation


##Moving to the QC
#first- compute nuclesome signal and store it back into the x.seurat object
x.seurat <- NucleosomeSignal(x.seurat)

#Second- calculating the  Transcriptional Start Site (TSS)
x.seurat <-TSSEnrichment(object = x.seurat, fast= FALSE)

#calculate the ratio of the fragments in the blaclisted region and the fraction of reads in the peak
x.seurat$blacklist_ratio <-x.seurat$blacklist_region_fragments/ x.seurat$peak_region_fragments
#calculate the fraction of reads in the peaks
x.seurat$pct_read_in_peaks <- x.seurat$peak_region_fragments/x.seurat$passed_filters*100

view(x.seurat)

#after viewing, just confirm if you have all the featured columns that you generated above
#Now let's move on to the visualization

#First get a list of the column names in the metadata
colnames(x.seurat@meta.data)

# 1) Let's make a scatter plot of two variables to see how they look
a1 <- DensityScatter(x.seurat, x = 'nCount_ATAC', y = 'TSS.enrichment',log_x = TRUE, quantiles = T)
#2)Let's make another density scatter plot for nucleosome signal vs TSS enrichment
a2<- DensityScatter(x.seurat, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = T)

#Let's visualize them together
CombinePlots(plots = list(a1, a2), ncol = 2)
#or just in a simple way 
a1 + a2
#Now let's make some violin plots to visualize the distribution of the metrics
VlnPlot(x.seurat, features = c('nCount_ATAC', 'nucleosome_signal', 'TSS.enrichment', 'pct_read_in_peaks', 'blacklist_ratio'), ncol = 5)

# After you view the plots, you can decide on the thresholds for filtering. Decide on the thresholds based on the distribution of the metrics
#Don't be too stringent with the thresholds as you can always filter more later if needed

# Filtering poor quality cells based on the metrics
x.seurat <- subset(
  x.seurat,
  subset = nCount_ATAC < 50000 &
    nCount_ATAC > 1000 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    pct_read_in_peaks > 15 &
    blacklist_ratio < 0.05
)

#After filtering, you can check the dimensions of the object to see how many cells you have left
dim(x.seurat)
#You can also visualize the metrics again to see how they look after filtering
VlnPlot(x.seurat, features = c('nCount_ATAC', 'nucleosome_signal', 'TSS.enrichment', 'pct_read_in_peaks', 'blacklist_ratio'), ncol = 5)

#Now let's move on to normalization and linear dimensional reduction
#First we will run the term frequency-inverse document frequency (TF-IDF) normalization
x.seurat <- RunTFIDF(x.seurat)

#Let's find the top features using the function FindTopFeatures
x.seurat <- FindTopFeatures(x.seurat, min.cutoff = 'q0')
#Now we will run Singular Value Decomposition (SVD) on the TF-IDF matrix (dimensiolity reduction)
x.seurat <- RunSVD(
  object = x.seurat,
  assay = 'ATAC',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

#The first componenet of the reduced dimension representation often captures the technical variation which is the sequencing depth
#rather than the biological variation. And this can be evaluated by visualizing the depth correlation 
DepthCor(x.seurat)
#You can see that the first component is highly correlated with the sequencing depth. So we will use the second component onwards for further analysis

#Now we will run UMAP for non-linear dimension reduction
x.seurat <- RunUMAP(object = x.seurat, reduction = 'lsi', dims = 2:30) #excluding the first component
#Now we will find the neighbors and clusters
x.seurat <- FindNeighbors(object = x.seurat, reduction = 'lsi', dims = 2:30)
x.seurat <- FindClusters(object = x.seurat, verbose = FALSE, algorithm = 3)
#Now let's visualize the clusters using UMAP
DimPlot(object = x.seurat, label = TRUE) + NoLegend()






