import os
import tacco as tc
import warnings

[ warnings.filterwarnings('ignore',category=c) for c in (FutureWarning, RuntimeWarning, UserWarning) ]

def _methods():
    return ['TACCO','SVM','SingleR','RCTD','NMFreg','Tangram*','NovoSpaRc*','WOT','NNLS','TACCO w/ cosine metric','TACCO w/ cosine metric, log-normalization','TACCO w/o platform normalization', 'TACCO w/o multicenter', 'TACCO w/o bisection',
'OT_BC','OTP_BC','OT_CS','OT_CSln','OT_BC_MC','OT_BC_BS','TACCO (cat)','TACCO w/ multicenter(5)','TACCO w/ multicenter(20)','TACCO w/ bisection(2,2)','TACCO w/ bisection(2,3)','TACCO w/ bisection(2,4)','TACCO w/ bisection(8,2)','TACCO w/ bisection(8,3)','TACCO w/ bisection(8,4)','TACCO w/ bisection(4,2)','TACCO w/ bisection(4,4)','TACCO w/ epsilon(0.05)','TACCO w/ epsilon(0.5)','TACCO w/ lambda(1.0)'
,'TACCO w/ lambda(0.01)',
]

def method_color(method):
    methods = _methods()
    colors = tc.pl.get_default_colors(methods)
    colors['TACCO (cat)'] = colors['TACCO']
    if method in methods:
        return colors[method]
    else:
        return '#777777'

def find_path(path_relative_to_base):
    # The notebook expects to be executed either in the workflow directory or in the repository root folder
    result_path = path_relative_to_base
    if not os.path.exists(result_path):
        result_path = f'../../{result_path}'
    result_path = os.path.abspath(result_path)

    if not os.path.exists(result_path):
        raise ValueError(f'The path {path_relative_to_base!r} cannot be resolved!')
    
    return result_path
