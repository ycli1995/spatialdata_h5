
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import spatialdata as sd
import spatialdata_plot as sdp  # noqa: F401
import squidpy as sq
from shapely.geometry import Point
from spatialdata.models import Image2DModel, ShapesModel, TableModel

adata = sq.datasets.visium_hne_adata()

genes = adata.var.nlargest(200, "n_cells_by_counts").index
cells = adata.obs.nlargest(100, "total_counts").index

adata = adata[cells, genes]
print(adata)

adata.write_h5ad("visium_mini.h5ad")

