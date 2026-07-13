import warnings

warnings.filterwarnings("ignore")

import json
import os
import scanpy as sc
import numpy as np
import sys

CT = sys.argv[1]

adata = sc.read_h5ad("data/Perturbseq/" + CT + ".h5ad")
adata.obs["mt_outlier"] = (
    adata.obs["mitopercent"] > 0.3
)

adata.obs.mt_outlier.value_counts()
##already filtered out
adata = adata[~adata.obs.mt_outlier].copy()

##important for cNMF. filter out cells with 0 expression
sc.pp.filter_cells(adata, min_genes=100) 
sc.pp.filter_cells(adata, min_counts=100) 
sc.pp.filter_genes(adata, min_cells=10)


h5ad_dir='data/Perturbseq/filtered_data/'

adata.write(h5ad_dir + CT + ".h5ad")
