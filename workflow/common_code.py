import os
import tacco as tc
import warnings

[ warnings.filterwarnings('ignore',category=c) for c in (FutureWarning, RuntimeWarning, UserWarning) ]

def _methods():
    return ['TACCO','SVM','SingleR','RCTD','NMFreg','Tangram*','NovoSpaRc*','WOT','NNLS','TACCO with cosine metric','TACCO with cosine metric, log-normalization','TACCO without platform normalization', 'TACCO without multi_center', 'TACCO without bisection', 'OT_BC','OTP_BC','OT_CS','OT_CSln','OT_BC_MC','OT_BC_BS','TACCO (cat)','TACCO with multi_center(5)','TACCO with multi_center(20)','TACCO with bisection(2,2)','TACCO with bisection(2,3)','TACCO with bisection(2,4)','TACCO with bisection(8,2)','TACCO with bisection(8,3)','TACCO with bisection(8,4)','TACCO with bisection(4,2)','TACCO with bisection(4,4)','TACCO with epsilon(0.05)','TACCO with epsilon(0.5)','TACCO with lambda(1.0)','TACCO with lambda(0.01)',
]

def method_color(method):
    methods = _methods()
    colors = tc.pl.get_default_colors(methods)
    colors['TACCO (cat)'] = colors['TACCO']
    if method in methods:
        return colors[method]
    else:
        return '#777777'

def method_style(method):
    if method in ['TACCO (cat)', 'SingleR', 'SVM',]:
        return 'dashed'
    else:
        return 'solid'

def find_path(path_relative_to_base, create_if_not_existent=False):
    workflow_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(workflow_dir)
    
    result_path = os.path.join(base_dir,path_relative_to_base)
    
    if not os.path.exists(result_path):
        if create_if_not_existent:
            os.makedirs(result_path)
        else:
            raise ValueError(f'The path {path_relative_to_base!r} cannot be resolved!')
    
    return result_path
