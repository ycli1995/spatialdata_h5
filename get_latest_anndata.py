
import scanpy as sc
import anndata as ad

adata = sc.read("pbmc_small.h5ad")

sc.pp.neighbors(adata, use_rep = "pca")
sc.tl.umap(adata, key_added = "umap")
sc.tl.leiden(adata, flavor="igraph", n_iterations=2, resolution = 0.1)

sc.pl.umap(adata, color=["leiden"], show = False)

sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")

print(adata)

adata.write_h5ad("pbmc_small.h5ad")

