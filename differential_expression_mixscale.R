# Importing libraries  
library(tidyverse)
library(readxl)
library(writexl)
library(qs2)
library(Seurat)
library(Mixscale)
library(DropletUtils)
library(SingleCellExperiment)
library(gprofiler2)


# Setting the file path to the 10X files
dir_10x <- "/mnt/mass_storage1/mass_storage_projects/compass/"



# Reading in the processed seurat objects
so <- qs_read(file = paste0(dir_10x, "seurat_objects/2_integrated_Mixscale_filt.qs2"))
meta <- so@meta.data



# Filtering based on a command line argument to split up the iNeuron dataset into smaller pieces
## Command line arguments should be parsed as follows: 
### 1. starting gRNA index
### 2. ending gRNA index
### 3. exported filename
#args <- commandArgs()
#start_index <- as.numeric(args[6])
#end_index <- as.numeric(args[7])
#gRNA_PRTBs_filt <- gRNA_PRTBs[start_index:end_index]
#so_small <- subset(so, NT_gRNA %in% c(gRNA_PRTBs_filt, "NonTargeting"))

# Clearing memory
#rm(so)
#gc()



# Getting a list of perturbations per gRNA
#gRNA_PRTBs <- sort(unique(so$NT_gRNA))
#gRNA_PRTBs <- gRNA_PRTBs[gRNA_PRTBs != "NonTargeting"]



# Running differential expression per gRNA with Mixscale
#de_res_gRNA <- Mixscale::Run_wmvRegDE(object = so,   #so or so_small
#                                      assay = "RNA", 
#                                      slot = "counts",
#                                      labels = "NT_gRNA", 
#                                      nt.class.name = "NonTargeting",
#                                      PRTB_list = gRNA_PRTBs,  #gRNA_PRTBs or gRNA_PRTBs_filt
#                                      logfc.threshold = 0,
#                                      split.by = "cell_type",
#                                      verbose = TRUE,
#                                      full.results = FALSE)


# Getting a list of perturbations per target gene
gene_PRTBs <- sort(unique(so$target_gene))
gene_PRTBs <- gene_PRTBs[gene_PRTBs != "NonTargeting"]


# Running differential expression per filtered gene with Mixscale
de_res_gene <- Mixscale::Run_wmvRegDE(object = so,   #so or so_filt
                                      assay = "RNA", 
                                      slot = "counts",
                                      labels = "target_gene", 
                                      nt.class.name = "NonTargeting",
                                      PRTB_list = gene_PRTBs,
                                      logfc.threshold = 0,
                                      split.by = "cell_type",
                                      verbose = TRUE,
                                      full.results = FALSE)

# Saving the weighted differential expression per target gene list
qs_save(object = de_res_gene, file = "integrative_tables/weighted_Mixscale_DEGs_per_target-gene.qs2")



# Saving the output
print("Done with differential expression - saving output DEG list!")
#export_filename <- as.character(args[8])
export_filename <- "integrated_Mixscale_DEGs"
qs_save(object = de_res_gRNA, file = paste0("integrative_tables/", export_filename, ".qs2"))



# Testing output
#de_res_gRNA <- qs_read("iNeuron_gRNA_DEGs/iNeuron_DEGs_list_505_1008.qs2")
#de <- de_res_gRNA[["KDM6A-g1"]]































