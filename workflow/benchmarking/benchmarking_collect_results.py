import pandas as pd
import os
import json

collection_csv = snakemake.output['collection_csv']

evaluation_inputs = [ fn for fn in snakemake.input if os.path.basename(fn).startswith('evaluation') ]
resources_inputs = [ fn for fn in snakemake.input if os.path.basename(fn).startswith('resources') ]

evaluation_dfs = []
for fn in evaluation_inputs:
    try:
        evaluation_dfs.append(pd.read_json(fn,orient='index',convert_dates=False).T)
    except:
        # Most likely the annotation failed due to insufficient resources - which is all we need to know for the benchmark
        pass
resources_dfs = []
for fn in resources_inputs:
    try:
        resources_dfs.append(pd.read_json(fn,orient='index',convert_dates=False).T)
    except:
        # Most likely the annotation failed due to insufficient resources - which is all we need to know for the benchmark
        pass

evaluation_df = pd.concat(evaluation_dfs)
resources_df = pd.concat(resources_dfs)

result_df = pd.merge(evaluation_df,resources_df)

# adjust the names of the methods to work around restrictions imposed by the file system etc.
result_df['method'] = result_df['method'].map({
    'TACCO': 'TACCO',
    'NMFreg': 'NMFreg',
    'RCTD': 'RCTD',
    'WOT': 'WOT',
    'Tangram_ast': 'Tangram*',
    'NovoSpaRc_ast': 'NovoSpaRc*',
    'TACCO_cat': 'TACCO (cat)',
    'SingleR': 'SingleR',
    'SVM': 'SVM',
    'TACCO_without_platform_normalization': 'TACCO without platform normalization',
    'TACCO_without_multicenter': 'TACCO without multi_center',
    'TACCO_without_bisection': 'TACCO without bisection',
    'TACCO_with_multicenter_5': 'TACCO with multi_center(5)',
    'TACCO_with_multicenter_20': 'TACCO with multi_center(20)',
    'TACCO_with_bisection_2_2': 'TACCO with bisection(2,2)',
    'TACCO_with_bisection_2_4': 'TACCO with bisection(2,4)',
    'TACCO_with_bisection_8_2': 'TACCO with bisection(8,2)',
    'TACCO_with_bisection_8_4': 'TACCO with bisection(8,4)',
    'TACCO_with_epsilon_0.05': 'TACCO with epsilon(0.05)',
    'TACCO_with_epsilon_0.5': 'TACCO with epsilon(0.5)',
    'TACCO_with_lambda_1.0': 'TACCO with lambda(1.0)',
    'TACCO_with_lambda_0.01': 'TACCO with lambda(0.01)',
})

os.makedirs(os.path.dirname(collection_csv), exist_ok=True)
result_df.to_csv(collection_csv,index=False)
