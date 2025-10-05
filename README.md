The preprocessing_QC file includes steps that start from loading the single cell ATAC seq data.
I used one data from 10X genomic. The link to the data is here https://www.10xgenomics.com/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-v-1-0-1-1-standard-1-2-0
Before beginning to analyse ATAC seq data, you should make sure that you have the peak by cell matrix data (if it is from 10x then it can be in the form of HDF5, metadata
(which is usually the per barcode metrics file) in the form of csv, and the fragments file (gz format).

Once you have all this you can use the script that is uploaded to run the general ATAC seq workflow pipeline. Additional comments on the steps can be found within the R script itself.
