import anndata as ad
import tacco as tc
import pandas as pd
import os
import json

import sample_data

dataset = snakemake.wildcards['dataset']
method = snakemake.wildcards['method']
Nref = snakemake.wildcards['Nref']
Ndat = snakemake.wildcards['Ndat']
Ncore = snakemake.wildcards['Ncore']
MB = snakemake.wildcards['MB']

evaluation_out = snakemake.output['evaluation_out']

try:
    annotation = pd.read_pickle(snakemake.input['annotation_in'])
    data = ad.read(snakemake.input['data_in'])
    with open(snakemake.input['keys_in']) as f:
        keys = json.load(f)
except:
    # Most likely the annotation failed due to insufficient resources - which is all we need to know for the benchmark
    # Therefore we just satisfy the output file requirements and go on, to be able to run the workflow further down
    open(evaluation_out, 'a').close()
    exit(0) # proclaim instant success

data = sample_data.sample_data(data, Ndat)

data.obsm[method] = annotation

corr = tc.ev.compute_err(data, method, keys['true_key'], err_method='corr')[method]
L1 = tc.ev.compute_err(data, method, keys['true_key'], err_method='lp', p=1)[method]
L2 = tc.ev.compute_err(data, method, keys['true_key'], err_method='lp', p=2)[method]
mc = tc.ev.compute_err(data, method, keys['true_key'], err_method='max_correct')[method]

evaluation_result = {}
evaluation_result['corr'] = corr
evaluation_result['L1'] = L1
evaluation_result['L2'] = L2
evaluation_result['mc'] = mc
evaluation_result['dataset'] = dataset
evaluation_result['method'] = method
evaluation_result['Nref'] = int(Nref)
evaluation_result['Ndat'] = int(Ndat)
evaluation_result['Ncore'] = int(Ncore)
evaluation_result['MB'] = int(MB)

os.makedirs(os.path.dirname(evaluation_out), exist_ok=True)
with open(evaluation_out, 'w') as f:
    json.dump(evaluation_result, f, indent=4)
