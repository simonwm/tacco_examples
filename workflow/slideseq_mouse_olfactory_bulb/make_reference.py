import urllib.request
import io
import gzip

import numpy as np
import scipy
import pandas as pd
import anndata as ad
import tacco as tc

def read_csv_from_gzipped_url(url,sep=','):
    with urllib.request.urlopen(url) as f:
        with gzip.GzipFile(mode='r', fileobj=f) as gz:
            with io.TextIOWrapper(gz) as gzt:
                table = pd.read_csv(gzt,sep=sep)
    return table

single_cell_expression_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121891/suppl/GSE121891_OB_6_runs.raw.dge.csv.gz'
single_cell_expression = read_csv_from_gzipped_url(single_cell_expression_url)

single_cell_metadata_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121891/suppl/GSE121891_OB_metaData_seurat.csv.gz'
single_cell_metadata = read_csv_from_gzipped_url(single_cell_metadata_url)

cluster2type = {
    'OEC': [ f'OEC{i+1}' for i in range(5)],
    'N': [ f'N{i+1}' for i in range(16)],
    'Astro': [ f'Astro{i+1}' for i in range(3)],
    'EC': [ f'EC{i+1}' for i in range(2)],
    'MicroG': [ f'MicroG{i+1}' for i in range(3)],
    'Mural': [ f'Mural{i+1}' for i in range(2)],
    'Mes': [ f'Mes{i+1}' for i in range(2)],
}
type2long = {
    'olfactory ensheathing cell-based (Sox10+)': 'OEC',
    'neuronal (Syt1+/Tubb3+)': 'N',
    'astrocytic (Gfap+)': 'Astro',
    'endothelial (Slco1c1+)': 'EC',
    'microglia (Aif1+/Siglech+)': 'MicroG',
    'myelinating-oligodendrocyte-based (Mag+)': 'MyOligo',
    'mural (Pdgfrb+)': 'Mural',
    'mesenchymal': 'Mes',
    'monocyte (Aif1+/Cd74+)': 'Mono',
    'macrophage (Aif1+/Cd52+)': 'Mφ',
    'oligodendrocyte-precursor-based (Olig2+)': 'OPC',
    'red blood cell (Aif1+/Hba-a1+)': 'RBCs',
}

tc.utils.merge_annotation(single_cell_metadata, 'ClusterName', cluster2type, 'type');
tc.utils.merge_annotation(single_cell_metadata, 'type', type2long, 'long');

single_cell_metadata = single_cell_metadata.set_index(single_cell_metadata.columns[0])
single_cell_expression = single_cell_expression[single_cell_metadata.index]

reference = ad.AnnData(single_cell_expression.T,obs=single_cell_metadata, dtype=np.int16)

reference.X = scipy.sparse.csr_matrix(reference.X)

reference.write(snakemake.output[0], compression='gzip', compression_opts=9)
