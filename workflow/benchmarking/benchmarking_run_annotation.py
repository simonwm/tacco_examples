annotation_out = snakemake.output['annotation_out']
resources_out = snakemake.output['resources_out']

attempt = snakemake.resources['attempt']

if attempt > 1: # we have been here before...
    # Most likely the last try failed due to insufficient resources - which is all we need to know for the benchmark
    # Therefore we just satisfy the output file requirements and go on, to be able to run the workflow further down
    open(annotation_out, 'a').close()
    open(resources_out, 'a').close()
    exit(0) # proclaim instant success

import scanpy as sc
import tacco as tc
import numpy as np
import pandas as pd
import os
import pickle
import h5py
import json

import sample_data

dataset = snakemake.wildcards['dataset']
method = snakemake.wildcards['method']
Nref = snakemake.wildcards['Nref']
Ndat = snakemake.wildcards['Ndat']
Ncore = snakemake.wildcards['Ncore']
MB = snakemake.wildcards['MB']

reference = sc.read(snakemake.input['reference_in'])
data = sc.read(snakemake.input['data_in'])
with open(snakemake.input['keys_in']) as f:
    keys = json.load(f)

reference = sample_data.sample_data(reference, Nref)
data = sample_data.sample_data(data, Ndat)

def get_method_config(method=None):

    method_kwargs = {
        'TACCO_without_bisection': {'method': 'OT', 'metric':'bc', 'multi_center': 10,'bisections':0,'platform_iterations': 0,},
        'TACCO_cat': {'method': 'OT', 'metric':'bc', 'multi_center': 10,'bisections':0,'platform_iterations': 0,'max_annotation':1,},
        'TACCO_without_multicenter': {'method': 'OT', 'metric':'bc', 'platform_iterations': 0,'bisections':4, 'bisection_divisor':3, },
        'TACCO_without_platform_normalization': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'bisections':4, 'bisection_divisor':3, 'platform_iterations': -1, },
        'TACCO': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'platform_iterations': 0,'bisections':4, 'bisection_divisor':3, },
        'TACCO_with_multicenter_5': {'method': 'OT', 'metric':'bc', 'multi_center': 5, 'platform_iterations': 0,'bisections':4, 'bisection_divisor':3, },
        'TACCO_with_multicenter_20': {'method': 'OT', 'metric':'bc', 'multi_center': 20, 'platform_iterations': 0,'bisections':4, 'bisection_divisor':3, },
        'TACCO_with_bisection_2_2': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'platform_iterations': 0,'bisections':2, 'bisection_divisor':2, },
        'TACCO_with_bisection_2_4': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'platform_iterations': 0,'bisections':2, 'bisection_divisor':4, },
        'TACCO_with_bisection_8_2': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'platform_iterations': 0,'bisections':8, 'bisection_divisor':2, },
        'TACCO_with_bisection_8_4': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'platform_iterations': 0,'bisections':8, 'bisection_divisor':4, },
        'TACCO_with_epsilon_0.05': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'platform_iterations': 0,'bisections':4, 'bisection_divisor':3,  'epsilon': 0.05},
        'TACCO_with_epsilon_0.5': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'platform_iterations': 0,'bisections':4, 'bisection_divisor':3,  'epsilon': 0.5},
        'TACCO_with_lambda_1.0': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'platform_iterations': 0,'bisections':4, 'bisection_divisor':3,  'lamb': 1.0},
        'TACCO_with_lambda_0.01': {'method': 'OT', 'metric':'bc', 'multi_center': 10, 'platform_iterations': 0,'bisections':4, 'bisection_divisor':3,  'lamb': 0.01},
        'NMFreg':{'method': 'NMFreg', 'K': 30, 'platform_iterations': -1, },
        'WOT':{'method': 'WOT', 'platform_iterations': -1, },
        'RCTD':{'method': 'RCTD', 'conda_env': os.path.abspath(snakemake.input['RCTD_env_link']), 'platform_iterations': -1, },
        'Tangram_ast':{'method':'tangram', 'conda_env': os.path.abspath(snakemake.input['tangram_env_link']), 'platform_iterations': -1, 'cluster_mode': False, },
        'NovoSpaRc_ast':{'method':'novosparc', 'platform_iterations': -1, 'n_hvg': 1000, },
        'SingleR':{'method':'SingleR', 'genes':'de', 'de_method':'classic', 'aggr_ref':False, 'conda_env': os.path.abspath(snakemake.input['SingleR_env_link']), 'platform_iterations': -1, },
        'SVM':{'method':'svm', 'platform_iterations': -1, },
    }

    if method in method_kwargs:
        return method_kwargs[method]
    else:
        raise ValueError(f'The method {method!r} is unknown!')

benchmark_result = tc.benchmarking.benchmark_annotate(data,reference,annotation_key=keys['ref_key'], **get_method_config(method), assume_valid_counts=True,)

os.makedirs(os.path.dirname(annotation_out), exist_ok=True)
benchmark_result['annotation'].to_pickle(annotation_out)

del benchmark_result['annotation']
benchmark_result['dataset'] = dataset
benchmark_result['method'] = method
benchmark_result['Nref'] = int(Nref)
benchmark_result['Ndat'] = int(Ndat)
benchmark_result['Ncore'] = int(Ncore)
benchmark_result['MB'] = int(MB)
os.makedirs(os.path.dirname(resources_out), exist_ok=True)
with open(resources_out, 'w') as f:
    json.dump(benchmark_result, f, indent=4)
