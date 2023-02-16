import anndata as ad
import numpy as np
import pandas as pd
import os
import time
import json

import tacco as tc

anno_csv = snakemake.input['anno_csv']
output_anno_csv = snakemake.output['seg_csv']
output_time_json = snakemake.output['time_json']

# load single molecule annotations
rna_coords = pd.read_csv(anno_csv)

# run single molecule segmentation
seg_col = f'seg'
timings = {}
for anno_col in ['tacco', 'baysor', 'ssam']:
    annseg_col = '%s_%s' % (anno_col, seg_col)
    print(f'segmenting single molecules: {anno_col}')
    start = time.time()
    tc.tl.segment(rna_coords, distance_scale=2.0, max_size=1600, result_key=annseg_col,
              position_scale=10.0,position_range=2,
              annotation_key=anno_col,annotation_distance=None,
                 )
    timings[annseg_col] = time.time() - start

# export the result
os.makedirs(os.path.dirname(output_anno_csv), exist_ok=True)
rna_coords.to_csv(output_anno_csv, index=False)
os.makedirs(os.path.dirname(output_time_json), exist_ok=True)
with open(output_time_json, 'w') as f:
    json.dump(timings,f)