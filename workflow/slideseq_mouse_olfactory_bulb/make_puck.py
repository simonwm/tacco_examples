import numpy as np
import scipy
import scipy.io
import pandas as pd
import anndata as ad

bead_locations = pd.read_csv(snakemake.input['location'])
bead_locations.columns = ['barcodes','x','y']

bead_expression = scipy.sparse.csr_matrix(scipy.io.mmread(snakemake.input['expression']))

gene_names = pd.read_csv(snakemake.input['gene_names'])
gene_names = gene_names[gene_names.columns[-1]].to_numpy()

puck = ad.AnnData(bead_expression.T,obs=bead_locations.set_index('barcodes'),var=pd.DataFrame(index=gene_names), dtype=np.int16)

puck = puck[puck.X.sum(axis=1).A.flatten() >= 50].copy()

puck.write(snakemake.output[0], compression='gzip', compression_opts=9)
