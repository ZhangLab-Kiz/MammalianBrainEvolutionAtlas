##########################################################################################
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad 
import os
import bbknn

results_file = 'singlecell.h5ad' 
merge_file = 'merge.h5ad'

base_dir = "."
ids = ['001', '002', '003']
subdirs = ['CC', 'CE', 'CS', 'TH', 'HY', 'Hi']

data_list = []

for id in ids:
    adata_list = []
    for subdir in subdirs:
        file_path = os.path.join(base_dir, id, subdir, 'filtered_cell_gene_matrix')
        adata = sc.read_10x_mtx(file_path, var_names='gene_symbols', cache=True)
        adata.obs['Organ'] = subdir  
        adata_list.append(adata)
        print(f"Loaded data from {file_path}")
    

    data = adata_list[0].concatenate(adata_list[1:])
    data.obs['ID'] = id  
    data_list.append(data)


adata_merged = data_list[0].concatenate(data_list[1:])
adata_merged.write(merge_file)
