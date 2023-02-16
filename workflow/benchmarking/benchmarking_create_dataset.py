import anndata as ad
import tacco as tc
import os
import json
import scipy.sparse

dataset = snakemake.wildcards['dataset']

reference_out = snakemake.output['reference_out']
data_out = snakemake.output['data_out']
keys_out = snakemake.output['keys_out']

if dataset == 'Mixture':
    
    os.symlink(os.path.abspath(snakemake.input['mouse_colon_scrnaseq']), os.path.abspath(reference_out)),
    
    reference = ad.read(snakemake.input['mouse_colon_scrnaseq'])
    
    capture_rate = 1.0
    bead_shape = 'gauss'
    ntdata_max = 10**4
    bead_size = 1.0

    ref_annotation_key = 'labels'
    tdata_annotation_key = 'reads_' + ref_annotation_key

    sdata = tc.tl.mix_in_silico(reference, type_key=ref_annotation_key, n_samples=ntdata_max, bead_shape=bead_shape, bead_size=bead_size, capture_rate=capture_rate,)
    sdata.obsm[tdata_annotation_key] /= sdata.obsm[tdata_annotation_key].to_numpy().sum(axis=1)[:,None]
    
    sdata.write(data_out, compression='gzip')
    
    keys = {
        'dataset': dataset,
        'ref_key': ref_annotation_key,
        'true_key': tdata_annotation_key,
    }
    with open(keys_out, 'w') as f:
        json.dump(keys, f, indent=4)
    
elif dataset == 'Dropout':
    from scsim import scsim   
 
    ngenes = 25000
    descale = 1.0
    ndoublets = 100
    K = 13
    nproggenes = 1000
    proggroups = [1,2,3,4]
    progcellfrac = .35
    ncells = 1500
    deprob = .025

    seed = 111

    deloc = 2.0

    # simulating true counts (in simulator.counts)
    simulator = scsim(ngenes=ngenes, ncells=ncells, ngroups=K, libloc=7.64, libscale=0.78,
                 mean_rate=7.68,mean_shape=0.34, expoutprob=0.00286,
                 expoutloc=6.15, expoutscale=0.49,
                 diffexpprob=deprob, diffexpdownprob=0., diffexploc=deloc, diffexpscale=descale,
                 bcv_dispersion=0.448, bcv_dof=22.087, ndoublets=ndoublets,
                 nproggenes=nproggenes, progdownprob=0., progdeloc=deloc,
                 progdescale=descale, progcellfrac=progcellfrac, proggoups=proggroups,
                 minprogusage=.1, maxprogusage=.7, seed=seed)
    simulator.simulate()

    drop_ref = ad.AnnData(scipy.sparse.csr_matrix(simulator.counts), obs=simulator.cellparams, var=simulator.geneparams)
    drop_ref.obs['group'] = drop_ref.obs['group'].astype('category')

    dropshape, dropmidpoint = simulator.fit_dropout()

    simulator.dropshape = dropshape
    simulator.dropmidpoint = -1.0
    simulator.simulate_dropouts()

    drop_data = ad.AnnData(scipy.sparse.csr_matrix(simulator.countswdrop), obs=simulator.cellparams, var=simulator.geneparams)
    drop_data.obs['group'] = drop_data.obs['group'].astype('category')

    drop_ref.write(reference_out, compression='gzip')
    
    drop_data.write(data_out, compression='gzip')
    
    keys = {
        'dataset': dataset,
        'ref_key': 'group',
        'true_key': 'group',
    }
    with open(keys_out, 'w') as f:
        json.dump(keys, f, indent=4)
    
elif dataset == 'Differentiation':
    
    os.symlink(os.path.abspath(snakemake.input['d4d6_h5ad']), os.path.abspath(reference_out)),
    
    os.symlink(os.path.abspath(snakemake.input['d2_h5ad']), os.path.abspath(data_out)),
    
    keys = {
        'dataset': dataset,
        'ref_key': 'Cell type annotation',
        'true_key': 'clone_fate',
    }
    with open(keys_out, 'w') as f:
        json.dump(keys, f, indent=4)
else:
    raise ValueError(f'The dataset {dataset!r} is not implemented!')
