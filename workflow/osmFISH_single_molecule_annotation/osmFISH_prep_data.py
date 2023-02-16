import scanpy as sc
import tacco as tc
import numpy as np
import pandas as pd
import os
import pickle
import tables

counts_loom = snakemake.input['counts_loom']
coords_hdf5 = snakemake.input['coords_hdf5']
seg_pkl = snakemake.input['seg_pkl']
output_counts_h5ad = snakemake.output['counts_h5ad']
output_coords_csv = snakemake.output['coords_csv']
output_coords_pxl_csv = snakemake.output['coords_pxl_csv']
output_centroid_profiles_csv = snakemake.output['centroid_profiles_csv']


# helper functions
def make_dataframe_from_coords_dic(coords_dic, label_key):
    coords = np.concatenate([v for v in coords_dic.values()],axis=0)
    keys = np.concatenate([np.full(len(v), k) for k,v in coords_dic.items()],axis=0)
    return pd.DataFrame({label_key:keys,'x':coords[:,0],'y':coords[:,1]})
def hybridization_annotation(hybridization_label):
    # annotation from Supplementary Table 2 from Codeluppi, S., Borm, L.E., Zeisel, A. et al. Spatial organization of the somatosensory cortex revealed by osmFISH. Nat Methods 15, 932–935 (2018). https://doi.org/10.1038/s41592-018-0175-z
    if hybridization_label.endswith('_Hybridization4'):
        return 'Partly_out_of_focus'
    if hybridization_label.startswith('Klk6_') or hybridization_label.startswith('Lum_'):
        return 'Low quality'
    if hybridization_label == 'Tbr1_Hybridization11':
        return 'Internal control'
    if hybridization_label.startswith('Cnr1') or hybridization_label.startswith('Plp1_') or hybridization_label.startswith('Vtn_'):
        return 'Repeat Round 4'
    return ''

# load segmented data
reference = sc.read_loom(counts_loom,sparse=False)
# get the scale
um_per_pixel = np.sqrt(reference.obs['size_um2']/reference.obs['size_pix']).mean()
reference.uns['um_per_pixel'] = um_per_pixel

# load segmentation data
with open(seg_pkl, 'rb') as f:
    seg_dic_all= pickle.load(f)
seg_df_all = make_dataframe_from_coords_dic(coords_dic=seg_dic_all, label_key='CellID')

# load single molecule data
with tables.open_file(coords_hdf5, mode="r") as h5file:
    channels = [ e.name for e in h5file if hasattr(e,'name') ]
    good_channels = [ c for c in channels if hybridization_annotation(c) not in ['Partly_out_of_focus','Low quality','Internal control'] ]
    
    #### match gene names
    sm_genes = pd.Series({ c: c.split('_')[0]  for c in good_channels })
    seg_genes = reference.var.index
    sm2seg_gene_matching = {sm:seg for sm,seg in zip(sorted(sm_genes), sorted(seg_genes))} # matching by sorting takes care of very simple cases only, e.g. renaming 'Kcnip'->'Kcnip2' and 'Tmem6'->'Tmem2' as needed here
    
    coords_dic = { sm2seg_gene_matching[sm_genes[c]]: h5file.root[c].read().astype(np.uint32) for c in good_channels }

rna_coords = make_dataframe_from_coords_dic(coords_dic=coords_dic, label_key='gene')

# create consistent hash for x,y coordinates
tc.utils.hash(seg_df_all,['x','y'],hash_key='hash',other=rna_coords, compress=False);

# get subset of segmentation assignments with RNA measurement
rna_coords_covered = rna_coords[rna_coords['hash'].isin(seg_df_all['hash'])]
seg_df_covered = seg_df_all[seg_df_all['hash'].isin(rna_coords_covered['hash'])]

# assign rna to segmented cell
rna_coords['cell'] = rna_coords['hash'].map(seg_df_covered.drop_duplicates('hash').set_index('hash')['CellID']) # drop_duplicates as the segmentation assigns some pixels/molecules to two cells
rna_coords.drop(columns=['hash'])
rna_coords['ClusterName'] = rna_coords['cell'].map(reference.obs['ClusterName']).fillna('')

# remove area used for stripping test
max_valid_x = rna_coords.query('ClusterName != "Excluded" & ClusterName != ""')['x'].max()
rna_coords = rna_coords.query(f'x <= {max_valid_x}')

# remove "Excluded" annotation
rna_coords.loc[rna_coords['ClusterName']=="Excluded",'ClusterName'] = ''
reference = reference[reference.obs['ClusterName'] != 'Excluded'].copy()

# optimize the output
rna_coords.loc[rna_coords['ClusterName']=='','ClusterName'] = np.nan
rna_coords['ClusterName'] = rna_coords['ClusterName'].astype('category')
rna_coords['cell'] = rna_coords['cell'].astype('category')
rna_coords['gene'] = rna_coords['gene'].astype(pd.CategoricalDtype(categories=reference.var.index, ordered=True))
reference.obs['ClusterName'] = reference.obs['ClusterName'].astype('category')

# export result in pixel coordinates and as csv for Baysor...
os.makedirs(os.path.dirname(output_coords_pxl_csv), exist_ok=True)
rna_coords.to_csv(output_coords_pxl_csv,index=False)

# transform to um
rna_coords[['x','y']] *= um_per_pixel
# reference coordinates are not needed and are therefore not transformed; it would also involve possible offsets and mirroring to get them aligned

# export the result
os.makedirs(os.path.dirname(output_counts_h5ad), exist_ok=True)
reference.write(output_counts_h5ad, compression='gzip')
os.makedirs(os.path.dirname(output_coords_csv), exist_ok=True)
rna_coords.to_csv(output_coords_csv,index=False)

# export cell type centroids as csv for Baysor...
one_hot = pd.get_dummies(reference.obs['ClusterName'])
one_hot.columns.name = 'ClusterName'
centroid_per_type = pd.DataFrame([np.mean(reference.X[which.to_numpy().astype(bool)], axis=0) for cname, which in one_hot.items()],index=one_hot.columns,columns=reference.var.index)
centroid_per_type.to_csv('temp.csv')
os.makedirs(os.path.dirname(output_centroid_profiles_csv), exist_ok=True)
centroid_per_type.to_csv(output_centroid_profiles_csv)
