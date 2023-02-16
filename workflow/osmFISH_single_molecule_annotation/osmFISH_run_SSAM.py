import anndata as ad
import numpy as np
import pandas as pd
import os
import time
import json

import ssam

counts_h5ad = snakemake.input['counts_h5ad']
coords_csv = snakemake.input['coords_csv']
output_anno_csv = snakemake.output['anno_csv']
output_time_json = snakemake.output['time_json']

# load segmented data
reference = ad.read(counts_h5ad)

# load single molecule data
rna_coords = pd.read_csv(coords_csv)

timings = {}
start = time.time()

# create SSAM data structure
width = rna_coords['x'].max() - rna_coords['x'].min() + 10 # add 5 um on every side to comfortably accomodate every molecule
height = rna_coords['y'].max() - rna_coords['y'].min() + 10
rna_coords_grouped = dict(iter(rna_coords.groupby('gene')))
genes = list(rna_coords_grouped.keys())
locations = [sub[['x','y']].to_numpy() for sub in rna_coords_grouped.values()]
ssamdata = ssam.SSAMDataset(genes, locations, width, height,)

# create SSAM analysis structure
ssamanalysis = ssam.SSAMAnalysis(ssamdata)

# run KDE
ssamanalysis.run_kde()
ssamanalysis.find_localmax(search_size=3, min_norm=0.027, min_expression=0.2, ) # use parameters from the ssam user guide
ssamanalysis.normalize_vectors_sctransform()

# calculate the cell type profiles
transformed = ssam.run_sctransform(reference[:,ssamdata.genes].X)[0].to_numpy()
one_hot = pd.get_dummies(reference[:,ssamdata.genes].obs['ClusterName'])
transformed_per_type = [np.mean(transformed[which.to_numpy().astype(bool)], axis=0) for cname, which in one_hot.items()]

# do the pixel-wise cell type mapping
ssamanalysis.map_celltypes(centroids=transformed_per_type)
pixel_type_map = ssamanalysis.dataset.celltype_maps

# map the pixel-wise cell type mapping back to single molecules
closest_pixel_coords = (rna_coords[['x','y']].to_numpy() + 0.5).astype(int).T
cell_type_id_per_molecule = pixel_type_map[closest_pixel_coords[0],closest_pixel_coords[1]].flatten()

cell_type_per_molecule = pd.Series(cell_type_id_per_molecule.astype(str))
cell_type_per_molecule.loc[cell_type_id_per_molecule < 0] = np.nan
cell_type_per_molecule.loc[cell_type_id_per_molecule >= 0] = one_hot.columns.astype(str)[cell_type_id_per_molecule[cell_type_id_per_molecule >= 0]]

timings['ssam'] = time.time() - start

# export the result
os.makedirs(os.path.dirname(output_anno_csv), exist_ok=True)
cell_type_per_molecule.to_csv(output_anno_csv, index=False, header=['type'])
os.makedirs(os.path.dirname(output_time_json), exist_ok=True)
with open(output_time_json, 'w') as f:
    json.dump(timings,f)