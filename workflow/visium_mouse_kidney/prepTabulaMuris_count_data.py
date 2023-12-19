#!/usr/bin/env python
# coding: utf-8

import os
import gzip
import tarfile
import zipfile
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

zip_archive = snakemake.input['droplet_zip']
metadata_file = snakemake.input['metadata']
output_h5ad = snakemake.output['scrnaseq_h5ad']

with zipfile.ZipFile(zip_archive, "r") as z:
    for file in z.namelist():
            if file.startswith("droplet/Kidney-"):
                z.extract(file, "resources/visium_mouse_kidney/TabulaMuris_count_data/")
adatas = {}
for sample in os.listdir("resources/visium_mouse_kidney/TabulaMuris_count_data/droplet/"):
    print(sample)
    adatas[sample] = sc.read_10x_mtx(os.path.join("resources/visium_mouse_kidney/TabulaMuris_count_data/droplet/", sample))
    adatas[sample].obs_names = adatas[sample].obs_names.map(lambda x: sample.split("-")[1] + "_" + x.split("-")[0])

adata = ad.concat(adatas.values())
adata = adata[(adata.obs_names!="10X_P4_6_TGCTACCTCCTTTACA"),:] #'Barcode (10X_P4_P6_)TGCTACCTCCTTTACA' is not in annotations_dropelets.csv and would result in Nan row

metadata = pd.read_csv(metadata_file).set_index("cell")
for column in metadata.columns:
    adata.obs[column] = metadata[column].astype("category")

# Filter celltypes that have 0 counts
adata.obs["cell_ontology_class"] = adata.obs["cell_ontology_class"].cat.remove_unused_categories()

os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)
adata.write(output_h5ad, compression='gzip')
