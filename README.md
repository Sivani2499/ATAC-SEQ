The preprocessing_QC file includes steps that start from loading the single cell ATAC seq data.
I used one data from 10X genomic. The link to the data is here https://www.10xgenomics.com/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-v-1-0-1-1-standard-1-2-0
Before beginning to analyse ATAC seq data, you should make sure that you have the peak by cell matrix data (if it is from 10x then it can be in the form of HDF5, metadata
(which is usually the per barcode metrics file) in the form of csv, and the fragments file (gz format).

Once you have all this you can use the script that is uploaded to run the general ATAC seq workflow pipeline. Additional comments on the steps can be found within the R script itself.

If you are new to ATAC seq data analysis, here is a table that tells you what each term means (some of the new functions in ATAC se)

| Term                                                   | Simple meaning                                                                                                                                                                                                                             |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **ChromatinAssay**                                     | A container in Seurat that stores all ATAC-seq information — the peaks (open DNA regions), fragment locations, and counts per cell. It’s the starting point for ATAC data analysis, similar to how “RNA Assay” holds gene expression data. |
| **Nucleosome signal**                                  | Tells how strongly DNA fragments show nucleosome spacing patterns. Low values mean more open chromatin (better-quality nuclei). High values suggest over-digested or low-quality cells.                                                    |
| **TSS enrichment**                                     | Measures how many reads fall near transcription start sites. High enrichment means active, well-captured cells with open promoters.                                                                                                        |
| **Blacklist ratio**                                    | Fraction of reads mapping to “blacklisted” genomic regions known to cause noise (repeats, artefacts). Low is good.                                                                                                                         |
| **FRiP (Fraction of Reads in Peaks)**                  | Shows what portion of all reads land inside true peaks (accessible DNA). High FRiP means cleaner data and better signal-to-noise.                                                                                                          |
| **TF-IDF (Term Frequency–Inverse Document Frequency)** | A normalization borrowed from text mining. It down-weights common peaks and up-weights unique ones, making cells more comparable.                                                                                                          |
| **SVD (Singular Value Decomposition)**                 | A mathematical method to reduce big matrices into a few main patterns. In ATAC analysis, it helps find major axes of variation.                                                                                                            |
| **LSI (Latent Semantic Indexing)**                     | The lower-dimensional representation that comes out of SVD + TF-IDF. It summarizes how similar cells are based on accessibility patterns, ready for UMAP and clustering.                                                                   |


An easy flowchart to keep in mind while starting ATAC seq analysis pre-processing:

        ┌────────────────────────────┐
        │   Input 10x ATAC files     │
        │ (.h5, fragments.tsv.gz,    │
        │  singlecell.csv)           │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ Read counts (.h5)          │
        │ + fragment file            │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ CreateChromatinAssay       │
        │  (store peaks + fragments) │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ CreateSeuratObject (ATAC)  │
        │  + add metadata            │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ Add gene annotations       │
        │  (EnsDb → UCSC "chr" fix)  │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ Quality Control (QC)       │
        │ • NucleosomeSignal         │
        │ • TSSEnrichment            │
        │ • FRiP & blacklist ratio   │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ Visualize metrics (plots)  │
        │ Decide filtering cutoffs   │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ Filter poor-quality cells  │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ Normalize with TF-IDF      │
        │ + FindTopFeatures          │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ Run SVD → get LSI          │
        │ (use LSI_2:30 components)  │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ Run UMAP (non-linear DR)   │
        │ FindNeighbors + Clusters   │
        └────────────┬───────────────┘
                     │
                     ▼
        ┌────────────────────────────┐
        │ Visualize clusters (UMAP)  │
        │ Save object (.rds)         │
        └────────────────────────────┘


🧬 Why do we create a ChromatinAssay?

Think of it as building the “ATAC compartment” inside the Seurat object.
Seurat can handle different data types (RNA, protein, ATAC, etc.), but each lives in its own assay slot.
When you call:
chrom_assay <- CreateChromatinAssay(counts = counts, fragments = "fragments.tsv.gz", ...)

you’re telling Seurat:
“Here are the open-chromatin counts and fragment positions — please package them neatly so downstream QC, TF-IDF, and LSI functions know where to look.”
Without this step, Seurat wouldn’t understand where your ATAC data lives or how to calculate accessibility metrics.

