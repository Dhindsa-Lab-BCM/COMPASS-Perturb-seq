import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from cnmf import cNMF, Preprocess
from scipy.io import mmread
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns



# Loading count matrix and metadata for iPSCs
ipsc_folder = "/mnt/mass_storage1/mass_storage_projects/compass/seurat_objects/iPSC_processed_sparse/"
counts = mmread(ipsc_folder + "2_integrated_Mixscale_filt_iPSC_counts.mtx").tocsr().T  #transpose to Cells x Genes
obs = pd.read_csv(ipsc_folder + "2_integrated_Mixscale_filt_iPSC_metadata.csv", index_col=0)
var = pd.read_csv(ipsc_folder + "2_integrated_Mixscale_filt_iPSC_genes.csv", index_col=1)

# Assemble and save
adata = sc.AnnData(X=counts, obs=obs, var=var)
#adata.write_h5ad(ipsc_folder + "cNMF_input_iPSC.h5ad")


# cNMF Analysis
cnmf_obj = cNMF(output_dir=ipsc_folder, name="iPSC_noBatchCorrection2")
k_range = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 19, 21, 23, 25, 27, 29, 30, 35, 40, 45, 50, 55, 60, 100]
cnmf_obj.prepare(counts_fn=ipsc_folder + 'cNMF_input_iPSC.h5ad', 
                 components=np.arange(5, 15), #testing K from 5 to 14
                 n_iter=20, 
                 seed=67, 
                 num_highvar_genes=2000)
cnmf_obj.factorize(worker_i=0, total_workers=1)
cnmf_obj.combine()
cnmf_obj.k_selection_plot()  #this generates a PNG in your output directory to help you pick K