import anndata as ad
import numpy as np
import pandas as pd
import os
import time
import json

import tacco as tc

counts_h5ad = snakemake.input['counts_h5ad']
coords_csv = snakemake.input['coords_csv']
output_anno_csv = snakemake.output['anno_csv']
output_time_json = snakemake.output['time_json']

# load segmented data
reference = ad.read(counts_h5ad)

# load single molecule data
rna_coords = pd.read_csv(coords_csv)

# run the single molecule annotation for different parameter choices
mods = {
    'base': { },
    'pla0': { 'platform_iterations': 0, },
    'bisWO': { 'bisections': 0, },
    'bis22': { 'bisections': 2, 'bisection_divisor': 2, },
    'bis24': { 'bisections': 2, 'bisection_divisor': 4, },
    'bis82': { 'bisections': 8, 'bisection_divisor': 2, },
    'bis84': { 'bisections': 8, 'bisection_divisor': 4, },
    'mul5': { 'multi_center': 5, },
    'mul10': { 'multi_center': 10, },
    'mul20': { 'multi_center': 20, },
    'eps05': { 'epsilon': 0.5, },
    'eps005': { 'epsilon': 0.05, },
    'lam1': { 'lamb': 1.0, },
    'lam001': { 'lamb': 0.01, },
    'bin20': { 'bin_size': 20, },
    'bin5': { 'bin_size': 5, },
    'bin20': { 'bin_size': 20, },
    'nsh2': { 'n_shifts': 2, },
    'nsh4': { 'n_shifts': 4, },
}
timings = {}
for mod_name, mod in mods.items():
    if mod_name == 'base':
        anno_mod_col = f'tacco'
    else:
        anno_mod_col = f'tacco_{mod_name}'
    annotation_parameters = { 'bin_size': 10, 'n_shifts': 3, 'bisections': 4, 'bisection_divisor': 3, 'platform_iterations': -1, }
    for k,v in mod.items():
        annotation_parameters[k] = v
    print(f'annotating single molecules with mod {mod_name}')
    start = time.time()
    tc.tl.annotate_single_molecules(rna_coords, reference=reference, method='OT', annotation_key='ClusterName', result_key=anno_mod_col, **annotation_parameters)
    timings[anno_mod_col] = time.time() - start

# export the result
os.makedirs(os.path.dirname(output_anno_csv), exist_ok=True)
rna_coords.to_csv(output_anno_csv,index=False)
os.makedirs(os.path.dirname(output_time_json), exist_ok=True)
with open(output_time_json, 'w') as f:
    json.dump(timings,f)