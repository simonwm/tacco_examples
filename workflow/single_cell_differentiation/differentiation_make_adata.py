import anndata as ad
import numpy as np
import scipy
import pandas as pd
import os

normed_counts = snakemake.input['normed_counts']
gene_names = snakemake.input['gene_names']
clone_matrix = snakemake.input['clone_matrix']
metadata = snakemake.input['metadata']

output_d4d6_h5ad = snakemake.output['d4d6_h5ad']
output_d2_h5ad = snakemake.output['d2_h5ad']

# load data
normed_counts_mat = scipy.io.mmread(normed_counts).tocsr()
genes = pd.read_csv(gene_names, sep='\t',header=None).to_numpy().flatten()
clone_mat = scipy.io.mmread(clone_matrix).tocsr()
meta_df = pd.read_csv(metadata, sep='\t')

# create full adata
adata = ad.AnnData(normed_counts_mat, obs=meta_df, var=pd.DataFrame(index=genes), dtype=np.float32)

# optimize dtypes
adata.obs['Library'] = adata.obs['Library'].astype('category')
adata.obs['Time point'] = adata.obs['Time point'].astype(int)
adata.obs['Starting population'] = adata.obs['Starting population'].astype('category')
adata.obs['Cell type annotation'] = adata.obs['Cell type annotation'].astype('category')
adata.obs['Well'] = adata.obs['Well'].astype(int)

# assign clone_id
adata.obs['clone_id'] = (clone_mat @ np.arange(1,1+clone_mat.shape[1])) - 1

# subset to clones appearing in day 2 and in day 4 or 6
time_dummies = pd.get_dummies(adata.obs['Time point'])
clone_in_day = pd.DataFrame(clone_mat.T @ time_dummies,columns=time_dummies.columns).astype(bool)
everywhere_clones = clone_in_day[2] & (clone_in_day[4] | clone_in_day[6])
everywhere_clone_ids = everywhere_clones[everywhere_clones].index
adata_everywhere = adata[adata.obs['clone_id'].isin(everywhere_clone_ids)]

# split into reference (day 4 and 6) and test (day 2) data
day4_day6 = adata_everywhere.obs['Time point']>2
adata_46 = adata_everywhere[day4_day6].copy()
adata_2 = adata_everywhere[~day4_day6].copy()

# determine clone fate (ground truth) for the test data from the differentiated data
clone_fate = adata_46.obs[['Cell type annotation', 'clone_id']].value_counts().unstack().fillna(0).T
adata_2.obsm['clone_fate'] = clone_fate.loc[adata_2.obs['clone_id']].set_index(adata_2.obs.index)
adata_2.obsm['clone_fate'] /= adata_2.obsm['clone_fate'].sum(axis=1).to_numpy()[:,None]

# use only "final" fates, i.e. all except 'Undifferentiated'
types = ['Neutrophil','Erythroid','Monocyte','Meg','Mast','Baso','Lymphoid','Eos','Ccr7_DC','pDC']
d4d6 = adata_46[adata_46.obs['Cell type annotation'].isin(types)]
d2 = adata_2[adata_2.obsm['clone_fate'][types].sum(axis=1) == 1]
d2.obsm['clone_fate'] = d2.obsm['clone_fate'][types]

# export the result
os.makedirs(os.path.dirname(output_d4d6_h5ad), exist_ok=True)
d4d6.write(output_d4d6_h5ad, compression='gzip')
os.makedirs(os.path.dirname(output_d2_h5ad), exist_ok=True)
d2.write(output_d2_h5ad, compression='gzip')
