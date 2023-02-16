import pandas as pd
import os

tacco_anno_csv = snakemake.input['tacco_csv']
baysor_anno_csv = snakemake.input['baysor_csv']
ssam_anno_csv = snakemake.input['ssam_csv']
output_anno_csv = snakemake.output['anno_csv']

# load single molecule annotations
rna_coords = pd.read_csv(tacco_anno_csv)
baysor = pd.read_csv(baysor_anno_csv)
rna_coords['baysor'] = baysor['clusterannotation']
rna_coords['baysor_beg'] = baysor['cell']
ssam = pd.read_csv(ssam_anno_csv)
rna_coords['ssam'] = ssam['type']

# export the result
os.makedirs(os.path.dirname(output_anno_csv), exist_ok=True)
rna_coords.to_csv(output_anno_csv, index=False)
