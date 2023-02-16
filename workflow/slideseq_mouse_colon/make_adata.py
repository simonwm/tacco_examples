import gzip
import scipy.io
import pandas as pd
import anndata as ad
import numpy as np
import tacco as tc

with gzip.open(snakemake.input['Slideseq_raw_mtx'], 'r') as f:
    X = scipy.io.mmread(f).T
cells = pd.read_csv(snakemake.input['Slideseq_raw_cells'], sep='\t', header=None).to_numpy().flatten()
genes = pd.read_csv(snakemake.input['Slideseq_raw_genes'], sep='\t', header=None).to_numpy().flatten()
slideseq_SCP_raw = ad.AnnData(X.tocsr(),obs=pd.DataFrame(index=cells),var=pd.DataFrame(index=genes),dtype=np.float32)

scRNAseq_SCP = ad.read(snakemake.input['scRNAseq_h5ad'])

singlecellportal_metadata = pd.read_csv(snakemake.input['singlecellportal_metadata'], sep='\t', skiprows=[1], index_col='NAME')

slideseq_SCP_raw.obs['x'] = singlecellportal_metadata['x_spatial']
slideseq_SCP_raw.obs['y'] = singlecellportal_metadata['y_spatial']
slideseq_SCP_raw.obs['sample'] = singlecellportal_metadata['puck']
slideseq_SCP_raw.obs['State'] = singlecellportal_metadata['disease__ontology_label']
slideseq_SCP_raw.obs['SampleID'] = singlecellportal_metadata['biosample_id']

slideseq = slideseq_SCP_raw
scrnaseq = scRNAseq_SCP

slideseq = slideseq[slideseq.obs['State'].isin(['normal']) & slideseq.obs['sample'].isin(['2020-09-14_Puck_200701_21'])]
scrnaseq = scrnaseq[scrnaseq.obs['State'].isin(['normal'])]

# make adata conform to the expected conventions
slideseq.obs[['x','y']] /= 0.65 # go back from µm to pixel units
slideseq.X = slideseq.X.tocsc()
slideseq = slideseq[:,tc.sum(slideseq.X,axis=0) != 0].copy()

slideseq.write(snakemake.output['slideseq_h5ad'], compression='gzip', compression_opts=9)
scrnaseq.write(snakemake.output['scrnaseq_h5ad'], compression='gzip', compression_opts=9)
