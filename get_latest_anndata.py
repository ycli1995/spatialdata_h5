
import scanpy as sc
import anndata as ad
import numpy as np

import random
import string

def random_len_str(max_len=10):
    length = random.randint(1, max_len)  # 随机长度 1~10
    return ''.join(random.choices(string.ascii_letters, k=length))

def const_len_str(length=8):
    return ''.join(random.choices(string.ascii_letters, k=length))

def random_str_arr(rows, cols, str_fun):
    arr = np.empty((rows, cols), dtype=object)
    for i in range(rows):
        for j in range(cols):
            arr[i, j] = str_fun()
    return arr


adata = sc.read("pbmc_small.h5ad")

sc.pp.neighbors(adata, use_rep = "pca")
sc.tl.umap(adata, key_added = "umap")
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution = 0.1)

sc.pl.umap(adata, color=["leiden"], show = False)

sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

adata.obsm['const_string_array'] = random_str_arr(len(adata.obs_names), 2, const_len_str)
adata.obsm['random_string_array'] = random_str_arr(len(adata.obs_names), 2, random_len_str)

adata.obsm['random_string_array_nan'] = adata.obsm['random_string_array']
adata.obsm['random_string_array_nan'][[1, 3, 5], [0, 1]] = None

adata.obsm['const_string_array_nan'] = adata.obsm['const_string_array']
adata.obsm['const_string_array_nan'][[1, 3, 5], [0, 1]] = None

adata.obsm['umap_nan'] = adata.obsm['umap']
adata.obsm['umap_nan'][[1, 3, 5], [0, 1]] = None

adata.obsm['bool'] = np.random.choice([True, False], size=(len(adata.obs_names), 2))
adata.obsm['bool_nan'] = adata.obsm['bool']
adata.obsm['bool_nan'][[1, 3, 5], [0, 1]] = None

adata.obs['integer_nan'] = adata.obs['nCount_RNA']
adata.obs['integer_nan'][[1, 3, 5]] = None

adata.obs['string'] = np.random.choice(['True', 'False'], size=len(adata.obs_names))
adata.obs['string_nan'] = adata.obs['string']
adata.obs['string_nan'][[1, 3, 5]] = None

adata.obs['bool'] = np.random.choice([True, False], size=len(adata.obs_names))
adata.obs['bool_nan'] = adata.obs['bool']
adata.obs['bool_nan'][[1, 3, 5]] = None

print(adata)

adata.write_h5ad("pbmc_small.h5ad")

