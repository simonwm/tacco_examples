#!/usr/bin/env python
# coding: utf-8
import os
import anndata as ad
import squidpy as sq

# Sample #1 Adult Mouse Kidney FFPE
visium_path = snakemake.input['visium_dir_ffpe']
output_h5ad = snakemake.output['visium_ffpe_h5ad']

adata = sq.read.visium(visium_path)
adata.obs['x'] = adata.obsm['spatial'][:,0]
adata.obs['y'] = -adata.obsm['spatial'][:,1]
adata.obs['slide'] = 'slide_1'
adata.obs['slide'] = adata.obs['slide'].astype("category") 

# Scale coordinates
scale = 304.75/65 #px/µm
adata.obs.loc[:,"x"]=adata.obs.loc[:,"x"].divide(scale) 
adata.obs.loc[:,"y"]=adata.obs.loc[:,"y"].divide(scale)

os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
adata.write(output_h5ad, compression='gzip')

# Sample #2 Adult Mouse Kidney Coronal v.1.10
visium_path = snakemake.input['visium_dir_coronal']
output_h5ad = snakemake.output['visium_coronal_h5ad']

adata = sq.read.visium(visium_path)
adata.obs['x'] = adata.obsm['spatial'][:,0]
adata.obs['y'] = -adata.obsm['spatial'][:,1]
adata.obs['slide'] = 'slide_2'
adata.obs['slide'] = adata.obs['slide'].astype("category") 

# Scale coordinates
scale = 144.51/65 #px/µm
adata.obs.loc[:,"x"]=adata.obs.loc[:,"x"].divide(scale) 
adata.obs.loc[:,"y"]=adata.obs.loc[:,"y"].divide(scale)

os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
adata.write(output_h5ad, compression='gzip')
