
import scanpy as sc
import anndata as ad
import numpy as np

adata = sc.read("example.h5ad")
adata.raw = adata

print(adata)

adata.write("latest_adata.h5ad")

