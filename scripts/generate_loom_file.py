# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from scipy import io
from scipy import sparse


scenic_path = "/kimdata/livinlrg/scAnalysis/BobDataComb/SCENIC/"
os.chdir(scenic_path)

object_path = "/kimdata/livinlrg/scAnalysis/BobDataComb/Objects/"

counts=io.mmread(object_path + 'WS290_elegans_expression_small.mtx')
counts=sparse.csr_matrix(counts)
barcodes=pd.read_csv(object_path + 'WS290_elegans_barcodes_small.csv')
genes=pd.read_csv(object_path + 'WS290_elegans_genes_small.csv')
adata=sc.AnnData(counts.T)
adata.obs_names=barcodes['Barcode'].values
adata.var_names=genes['Gene'].values

scratch_path = "/scratch/livinlrg/"
os.chdir(scratch_path)

f_loom_path_scenic = scratch_path + "WS290_elegans_expression.loom"

# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

# os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
lp.create("WS290_elegans_expression.loom", adata.X.transpose(), row_attrs, col_attrs)

counts=io.mmread(object_path + 'WS290_briggsae_expression.mtx')
counts=sparse.csr_matrix(counts)
barcodes=pd.read_csv(object_path + 'WS290_briggsae_barcodes.csv')
genes=pd.read_csv(object_path + 'WS290_briggsae_genes.csv')
adata=sc.AnnData(counts.T)
adata.obs_names=barcodes['Barcode'].values
adata.var_names=genes['Gene'].values

scratch_path = "/scratch/livinlrg/"
os.chdir(scratch_path)

f_loom_path_scenic = scratch_path + "WS290_briggsae_expression.loom"

# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

# os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
lp.create("WS290_briggsae_expression.loom", adata.X.transpose(), row_attrs, col_attrs)

