import numpy as np
import pandas as pd
import tacco as tc
import scanpy as sc
import warnings
warnings.filterwarnings('ignore')
import matplotlib
matplotlib.rcParams["axes.titlelocation"] = 'left'
matplotlib.rcParams["axes.titlesize"] = 'small'
import scipy

# Initialize data variables
slideseq = None
singlecell = None
means = None
types = None
dirtcell = None
txO = None
tyO = None

# Initialize plotting parameters
colors_empty_celltype = None
colors_celltype = None
colors_0 = None
colors_1 = None
colors_2 = None
colors_region = None
scale = None
radius_celltype = None
cent = None
rad = None
txmin = None
tymin = None
txmax = None
tymax = None

def generate_data():
    lattice_spacing = 1
    nx = 5
    x = np.linspace(0, lattice_spacing * nx, nx)
    y = np.linspace(0, lattice_spacing * nx, nx)
    xv, yv = np.meshgrid(x, y)
    xv += y[:,None]/2
    yv *= 1 / np.sqrt(1.25)

    center_x = (xv[-1,-1] - xv[0,0])/2
    center_y = (yv[-1,-1] - yv[0,0])/2
    radius = np.sqrt((xv[:,0] - center_x)**2 + (yv[:,0] - center_y)**2).min() + lattice_spacing / 2

    xv,yv = xv.flatten(),yv.flatten()

    np.random.seed(42)
    xv += np.random.normal(0.0,lattice_spacing*0.15,xv.shape)
    yv += np.random.normal(0.0,lattice_spacing*0.15,yv.shape)

    within = np.sqrt((xv - center_x)**2 + (yv - center_y)**2) <= radius
    xv = xv[within]
    yv = yv[within]

    np.random.seed(42)
    p_0_width = 1 * lattice_spacing
    p_0_center_i = np.random.choice(len(xv))
    p_0_center_x = xv[p_0_center_i]
    p_0_center_y = yv[p_0_center_i]
    p_0 = np.exp(-((xv - p_0_center_x)**2 + (yv - p_0_center_y)**2) / (2 * p_0_width**2))
    p_1_width = 2 * lattice_spacing
    p_1 = (1 - p_0) * np.exp(-((xv - p_0_center_x)**2 + (yv - p_0_center_y)**2) / (2 * p_1_width**2))
    p_2 = (1 - p_0 - p_1) * 1
    p_0_flat = 0.0
    p_1_flat = 0.0
    p_2_flat = 0.0
    p_flat = p_0_flat + p_1_flat + p_2_flat
    p_0 *= 1 - p_flat
    p_1 *= 1 - p_flat
    p_2 *= 1 - p_flat
    p_0 += p_0_flat
    p_1 += p_1_flat
    p_2 += p_2_flat

    global slideseq
    slideseq = pd.DataFrame({'x':xv, 'y':yv})
    evidence = 20
    np.random.seed(42)
    slideseq[['f0','f1','f2']] = np.array([np.random.dirichlet((p0*evidence+1,p1*evidence+1,p2*evidence+1)) for p0,p1,p2 in zip(p_0,p_1,p_2) ])
    
    slideseq['tx'] = 0.5 * (2 * slideseq['f1'] + slideseq['f2']) / (slideseq['f0'] + slideseq['f1'] + slideseq['f2'])
    slideseq['ty'] = np.sqrt(3) * 0.5 * (slideseq['f2']) / (slideseq['f0'] + slideseq['f1'] + slideseq['f2'])
    
    type_totals = slideseq[['f0','f1','f2']].sum(axis=0).round().astype(int)
    n_cells = type_totals.sum()
    global singlecell, types
    singlecell = pd.DataFrame({'x':np.full(shape=n_cells,fill_value=np.nan,),'y':np.full(shape=n_cells,fill_value=np.nan,)})
    types = []
    for t,tn in type_totals.items():
        for i in range(tn):
            types.append(t)
    types = np.array(types)
    for t in type_totals.index:
        singlecell[t] = (types == t) * 1.0
    singlecell['tx'] = 0.5 * (2 * singlecell['f1'] + singlecell['f2']) / (singlecell['f0'] + singlecell['f1'] + singlecell['f2'])
    singlecell['ty'] = np.sqrt(3) * 0.5 * (singlecell['f2']) / (singlecell['f0'] + singlecell['f1'] + singlecell['f2'])

    np.random.seed(42)
    singlecell['tx'] += np.random.normal(0.0,0.1,singlecell.shape[0])
    singlecell['ty'] += np.random.normal(0.0,0.1,singlecell.shape[0])

    singlecell['type'] = types
    global means
    means = singlecell.groupby(types).mean()
    
    tx,ty = np.concatenate([slideseq['tx'],singlecell['tx'],means['tx']]),np.concatenate([slideseq['ty'],singlecell['ty'],means['ty']])
    global txO,tyO
    txO,tyO = increase_distances(tx,ty)
    slideseq['txO'],slideseq['tyO'] = txO[:len(slideseq)],tyO[:len(slideseq)]
    singlecell['txO'],singlecell['tyO'] = txO[len(slideseq):-len(means)],tyO[len(slideseq):-len(means)]
    means['txO'],means['tyO'] = txO[-len(means):],tyO[-len(means):]

    for im in range(len(means)):
        slideseq[f'tx{im}'],slideseq[f'ty{im}'] = shrink_distances(slideseq['tx'],slideseq['ty'],means.iloc[im]['tx'],means.iloc[im]['ty'])
    
    global dirtcell
    dirtcell = singlecell.copy()
    dirtcell['dirt'] = ''
    np.random.seed(42)
    for group in means.index:
        sub_frac = type_totals[type_totals.index != group]
        sub_frac /= sub_frac.sum()
        sel = dirtcell['type'] == group
        dirtcell.loc[sel,'dirt'] = np.random.choice(sub_frac.index, size=sel.sum(), replace=True, p=sub_frac.to_numpy())
    dirtweight = 0.40
    dirtcell['tx'] = dirtcell['tx'] * (1-dirtweight) + dirtcell['dirt'].map(lambda x: means.loc[x,'tx']) * dirtweight
    dirtcell['ty'] = dirtcell['ty'] * (1-dirtweight) + dirtcell['dirt'].map(lambda x: means.loc[x,'ty']) * dirtweight
    dirtcell['txO'],dirtcell['tyO'] = increase_distances(dirtcell['tx'],dirtcell['ty'])

def initialize_plotting_params():
    global colors_empty_celltype, colors_celltype, colors_0, colors_1, colors_2, colors_region, scale, radius_celltype, cent, rad, txmin, tymin, txmax, tymax 
    colors_empty_celltype = ['#aaaaaa']
    colors_celltype = ['#0f0','#f00','#00f',]
    colors_0 = [ c if i == 0 else '#fff0' for i,c in enumerate(colors_celltype) ]
    colors_1 = [ c if i == 1 else '#fff0' for i,c in enumerate(colors_celltype) ]
    colors_2 = [ c if i == 2 else '#fff0' for i,c in enumerate(colors_celltype) ]
    colors_region = {'0':'#66a','1':'#6a6','2':'#a66',}

    scale = 1.0
    radius_celltype = 0.4
    cent = np.array([(slideseq['x'].max()+slideseq['x'].min())/2, (slideseq['y'].max()+slideseq['y'].min())/2])
    rad = max(np.sqrt(((slideseq[['x','y']] - cent)**2).sum(axis=1))) + 1.1 * radius_celltype
    
    txmin = min(slideseq['tx0'].min(), slideseq['tx1'].min(), slideseq['tx2'].min(), txO.min())
    tymin = min(slideseq['ty0'].min(), slideseq['ty1'].min(), slideseq['ty2'].min(), tyO.min())
    txmax = max(slideseq['tx0'].max(), slideseq['tx1'].max(), slideseq['tx2'].max(), txO.max())
    tymax = max(slideseq['ty0'].max(), slideseq['ty1'].max(), slideseq['ty2'].max(), tyO.max())
        
def increase_distances(x,y):
    from scipy.optimize import minimize
    def potential(xy,xy0,nbig):
        x,y = xy[:len(xy)//2],xy[len(xy)//2:]
        x0,y0 = xy0[:len(xy0)//2],xy0[len(xy0)//2:]
        pot = ((x-x0)**2 + (y-y0)**2).sum() # use a quadratic potential to keep points where they are ...
        pot += 3 * ((x[-nbig:]-x0[-nbig:])**2 + (y[-nbig:]-y0[-nbig:])**2).sum() # ... especially "big" points
        dists = tc.utils.cdist(np.array([x,y]).T, metric='euclidean') # and increase the distances between the points ...
        scaling = np.ones(len(dists))
        scaling[-nbig:] *= 0.25 # ... especially for big points ...
        dists *= scaling[None,:]
        dists = dists[np.triu_indices(len(x),1)]
        pot += 0.0003 * (dists**(-1)).sum() # ... using a 1/distance potential
        return pot
    xy0 = np.concatenate([x,y],)
    res = minimize(potential, xy0, method='BFGS', args=(xy0,3)) #options={'disp': True},
    return res.x[:len(res.x)//2],res.x[len(res.x)//2:]

def shrink_distances(x,y,xs,ys):
    from scipy.optimize import minimize
    def potential(xy,xy0,xs,ys):
        x,y = xy[:len(xy)//2],xy[len(xy)//2:]
        x0,y0 = xy0[:len(xy0)//2],xy0[len(xy0)//2:]
        pot = 0.0 * ((x-x0)**2 + (y-y0)**2).sum() # use a weak quadratic potential to keep the memory of the origin of the points
        pot += ((x-xs)**2 + (y-ys)**2).sum() # use a strong quadratic potential to attract points to the attractor
        dists = tc.utils.cdist(np.array([x,y]).T, metric='euclidean') # and increase the distances between the points ...
        dists = dists[np.triu_indices(len(x),1)]
        pot += 0.0012 * (dists**(-1)).sum() # ... using a 1/distance potential
        dists = tc.utils.cdist(np.array([x,y]).T, np.array([[xs],[ys]]).T, metric='euclidean') # ... also for the attractor
        dists = dists.flatten()
        pot += 0.0024 * (dists**(-1)).sum() # ... using a 1/distance potential
        return pot
    xy0 = np.concatenate([x,y],)
    res = minimize(potential, xy0, method='BFGS', args=(xy0,xs,ys)) #options={'disp': True},
    return res.x[:len(res.x)//2],res.x[len(res.x)//2:]

def plot_slideseq(ax):
    for row in slideseq.itertuples():
        ax.pie([1.0],colors=colors_empty_celltype,radius=radius_celltype*0.1,center=(row.txO*scale,row.tyO*scale), frame=True)
    set_ax_options(ax)

def plot_scRef(ax):
    ax.scatter(singlecell['txO'], singlecell['tyO'],marker='^',c=pd.Series(types).astype('category').cat.codes.map(lambda x: colors_celltype[x]))
    set_ax_options(ax, title='Categorical annotation \n(e.g., single cell data)')

def plot_scRef_profiles(ax):
    for irow,row in enumerate(means.itertuples()):
        patch = matplotlib.patches.RegularPolygon((row.txO*scale,row.tyO*scale), 3, radius_celltype*0.2, ec='#0006', fc=colors_celltype[irow], clip_on=False)
        ax.add_patch(patch)

def plot_slideseq_frac(ax):
    for row in slideseq.itertuples():
        ax.pie(row[3:6],colors=colors_celltype,radius=radius_celltype*0.1,center=(row.txO*scale,row.tyO*scale), frame=True, wedgeprops = {'ec':'#0002'})
        
def plot_ambiguous_anno(ax):
    cats = pd.Series(dirtcell['type']).astype('category').dtype
    ax.scatter(dirtcell['txO'],dirtcell['tyO'],marker='^', c=pd.Series(dirtcell['type']).astype(cats).cat.codes.map(lambda x: colors_celltype[x]),
            edgecolors=pd.Series(dirtcell['dirt']).astype(cats).cat.codes.map(lambda x: colors_celltype[x]),)
    
def plot_continuous_anno(ax):
    base_colors = matplotlib.colors.to_rgba_array(colors_celltype)
    mixed_colors = tc.pl.mix_base_colors(slideseq[['f0','f1','f2']], base_colors).to_numpy()
    for row_color,row in zip(mixed_colors,slideseq.itertuples()):
        ax.pie([1.0],colors=[row_color],radius=radius_celltype*0.1,center=(row.txO*scale,row.tyO*scale), frame=True, wedgeprops = {'ec':'#0002'})
        
def plot_OT_anno(ax):
    dists = tc.utils.cdist(np.array([slideseq['txO'],slideseq['tyO']]).T, np.array([means['txO'],means['tyO']]).T, metric='euclidean')
    far_enough = np.arange(len(dists))[(dists > 0.3).all(axis=1)]
    selected = []
    for i in range(means.shape[0]):
        indices = np.argsort(dists[:,i])
        selected.append(indices[np.isin(indices,far_enough)][0])

    for i,row in enumerate(slideseq.iloc[selected].itertuples()):
        ax.pie(row[3:6],colors=colors_celltype,radius=radius_celltype*0.1,center=(row.txO*scale,row.tyO*scale), frame=True, wedgeprops = {'ec':'#0002'})

        for mi,mrow in enumerate(means.itertuples()):

            patch = matplotlib.patches.RegularPolygon((mrow.txO*scale,mrow.tyO*scale), 3, radius_celltype*0.2,  ec='#0006', fc=colors_celltype[mi], clip_on=False)
            ax.add_patch(patch)

            arrowprops = dict(arrowstyle="->",connectionstyle="arc", facecolor=colors_celltype[mi], edgecolor=colors_celltype[mi],
                              shrinkA=3, shrinkB=6, patchA=patch)

            ax.annotate('', xy=(row.txO*scale, row.tyO*scale), xycoords='data', ha="left", va="center", fontsize='small',
                    xytext=(mrow.txO*scale, mrow.tyO*scale), textcoords='data', arrowprops=arrowprops,)

def plot_object_split_total(ax):
    for row in slideseq.itertuples():
        ax.pie(row[3:6],colors=colors_0,radius=radius_celltype*0.1,center=(row.tx0*scale,row.ty0*scale), frame=True, wedgeprops = {'ec':'#0002'})
        ax.pie(row[3:6],colors=colors_1,radius=radius_celltype*0.1,center=(row.tx1*scale,row.ty1*scale), frame=True, wedgeprops = {'ec':'#0002'})
        ax.pie(row[3:6],colors=colors_2,radius=radius_celltype*0.1,center=(row.tx2*scale,row.ty2*scale), frame=True, wedgeprops = {'ec':'#0002'})

def plot_compositional_annotation(ax):
    plot_scRef_profiles(ax)
    plot_slideseq_frac(ax)
    set_ax_options(ax, title='Compositional annotation\n(e.g., Slide-Seq)')
    
def plot_full_annotated(ax):
    plot_scRef_profiles(ax)
    plot_slideseq_frac(ax)
    set_ax_options(ax, title='Full annotated dataset\n ')
    
def plot_ambiguous_categories(ax):
    plot_ambiguous_anno(ax)
    plot_scRef(ax)
    set_ax_options(ax, title='Ambiguous categories\n(e.g., ambient, droput)')
    
def plot_continuous_spectrum(ax):
    plot_scRef_profiles(ax)
    plot_continuous_anno(ax)
    set_ax_options(ax, title='Continuous spectrum\n(e.g., differentiation)')
    
def plot_common_space(ax):
    plot_scRef(ax)
    plot_slideseq(ax)
    set_ax_options(ax, title='Common space for all data\n(after platform normalization)')
    
def plot_reference_profiles(ax):
    plot_scRef(ax)
    plot_slideseq(ax)
    plot_scRef_profiles(ax)
    set_ax_options(ax, title='Reference profiles per class\n(multiple for subclustering)')
    
def plot_optimal_transport(ax):
    plot_slideseq(ax)
    plot_OT_anno(ax)
    set_ax_options(ax, title='Optimal transport of annotations\n(multiple rounds for bisectioning)')
    
def plot_object_splitting(ax):
    plot_scRef_profiles(ax)
    plot_object_split_total(ax)
    set_ax_options(ax, title='Annotations after object splitting in expression space\n(e.g., for categorical downstream analyses)')
    
def set_ax_options(ax, title=None):
    ax.set_aspect('equal')
    
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([])

    ax.text(0.38, -0.04, 'gene A', transform=ax.transAxes,)
    ax.text(-0.04, 0.39, 'gene B', transform=ax.transAxes,rotation='vertical',)
    
    if title is not None:
        ax.set_title(title)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    arrow_args = {'facecolor':'#000', 'head_width':0.1, 'overhang':0.3, 'length_includes_head':True, 'clip_on':False}
    ax.arrow(txmin-0.1, tymin-0.1, txmax-txmin+0.2, 0, **arrow_args) 
    ax.arrow(txmin-0.1, tymin-0.1, 0, tymax-tymin+0.2, **arrow_args)

def generate_adata():
    adata = sc.AnnData(np.ones(len(slideseq))[:,None], obs=slideseq,)
    adata.obsm['fractions'] = adata.obs[['f0', 'f1', 'f2']]
    return adata

def set_pl_co_occ_params(fig, adata):
    ax = fig.get_axes()[0]
    ax.set_title('Longe-range structure\ne.g., co-occurence\n')
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_xlabel(None)
    ax.text(0.30, -0.25, 'distance', transform=ax.transAxes,)
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_ylabel(None)
    ax.text(-0.20, 0.17, 'p(*|   ,distance)', transform=ax.transAxes,rotation='vertical',)
    ax.text(-0.20, 0.17, '      \u2588', transform=ax.transAxes,rotation='vertical',color=colors_celltype[0])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    xmin, xmax = ax.get_xlim() 
    ymin, ymax = ax.get_ylim() 

    txmin = adata.uns['fractions']['interval'][0]
    tymin = adata.uns['fractions']['composition'].min()
    txmax = adata.uns['fractions']['interval'][-1]
    tymax = adata.uns['fractions']['composition'].max()

    arrow_args = {'facecolor':'#000', 'overhang':0.3, 'length_includes_head':True, 'clip_on':False}
    ax.arrow(txmin-0.1, tymin-0.1, txmax-txmin+0.2, 0, **arrow_args, head_width=0.05*(ymax-ymin), head_length=0.075*(xmax-xmin), ) 
    ax.arrow(txmin-0.1, tymin-0.1, 0, tymax-tymin+0.2, **arrow_args, head_width=0.05*(xmax-xmin), head_length=0.075*(ymax-ymin), )

    ax.set_xlim(xmin-0.01, xmax+0.01)
    ax.set_ylim(ymin-0.01, ymax+0.01)
    
def set_pl_co_occ_matrix(fig):
    ax = fig.get_axes()[0]
    ax.set_title('Short-range structure\n(e.g., neighbourship z-score)\n')
    ax.set_xticklabels(['      ','      ','      '])
    ax.set_xlabel('center')
    ax.set_yticklabels(['    ','    ','    '])
    ax.set_ylabel('annotation')
    ax.text(-0.20, 0.14+0/3, '\u2588', transform=ax.transAxes,color=colors_celltype[0], fontsize='x-large')
    ax.text(-0.20, 0.14+1/3, '\u2588', transform=ax.transAxes,color=colors_celltype[1], fontsize='x-large')
    ax.text(-0.20, 0.14+2/3, '\u2588', transform=ax.transAxes,color=colors_celltype[2], fontsize='x-large')
    ax.text(0.10+0/3, -0.20, '\u2588', transform=ax.transAxes,color=colors_celltype[0], fontsize='x-large')
    ax.text(0.10+1/3, -0.20, '\u2588', transform=ax.transAxes,color=colors_celltype[1], fontsize='x-large')
    ax.text(0.10+2/3, -0.20, '\u2588', transform=ax.transAxes,color=colors_celltype[2], fontsize='x-large')
    
def plot_unannotated_data(ax, adata):
    title='Unannotated data \n(e.g., Slide-Seq)'
    set_ax_options_region(ax, adata, title)
    for row in slideseq.itertuples():
        ax.pie([1.0],colors=colors_empty_celltype,radius=radius_celltype*scale,center=(row.x*scale,row.y*scale), frame=True, wedgeprops = {'ec':'#0002'})
    
def plot_annotated_data(ax, adata):
    set_ax_options_region(ax, adata, title='Annotated data \n ')
    for row in slideseq.itertuples():
        ax.pie(row[3:6],colors=colors_celltype,radius=radius_celltype*scale,center=(row.x*scale,row.y*scale), frame=True, wedgeprops = {'ec':'#0002'})
    
def plot_object_split_green(ax, adata):
    set_ax_options_region(ax, adata, title='Annotation 1 in position space\nafter object splitting') 
    for row in slideseq.itertuples():
         ax.pie(row[3:6],colors=colors_0,radius=radius_celltype*scale,center=(row.x*scale,row.y*scale), frame=True, wedgeprops = {'ec':'#0002'})
    
def plot_object_split_red(ax, adata):
    set_ax_options_region(ax, adata, title='Annotation 2 in position space\nafter object splitting')
    for row in slideseq.itertuples():
        ax.pie(row[3:6],colors=colors_1,radius=radius_celltype*scale,center=(row.x*scale,row.y*scale), frame=True, wedgeprops = {'ec':'#0002'})
        
def plot_object_split_blue(ax, adata):
    set_ax_options_region(ax, adata, title='Annotation 3 in position space\nafter object splitting')
    for row in slideseq.itertuples():
        ax.pie(row[3:6],colors=colors_2,radius=radius_celltype*scale,center=(row.x*scale,row.y*scale), frame=True, wedgeprops = {'ec':'#0002'})
    
def plot_region_voronoi(ax, adata):
    set_ax_options_region(ax, adata)
    for row in slideseq.itertuples():
        ax.pie(row[3:6],colors=colors_celltype,radius=radius_celltype*scale,center=(row.x*scale,row.y*scale), frame=True, wedgeprops = {'ec':'#0002'})
    set_ax_options_region_voronoi(ax, adata)
    
def set_ax_options_region(ax, adata, title=None):
    if title is not None:
        ax.set_title(title)
    
    ax.add_patch(matplotlib.patches.Circle(cent, rad, ec='black', fc='white', lw=1, clip_on=False, zorder=0.8))
    
    ax.scatter(slideseq['x'],slideseq['y'], zorder=0.7)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_xlabel('spatial x')
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.set_ylabel('spatial y')
    
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
def set_ax_options_region_voronoi(ax, adata):
    ax.set_title('Regions from spatial and \nannotation information')
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()

    points = slideseq[['x','y']].to_numpy()
    points = np.append(points, [[0,1e3], [0,-1e3], [1e3,0], [-1e3,0]], axis=0)
    vor = scipy.spatial.Voronoi(points)
    for pi,ri in enumerate(vor.point_region):
        region = vor.regions[ri]
        if not -1 in region:
            polygon = [vor.vertices[i] for i in region]
            ax.fill(*zip(*polygon),c=colors_region[adata.obs.iloc[pi]['region']], zorder=0.9)
    
    outer = matplotlib.path.Path.circle(cent,1*rad)
    inner = matplotlib.path.Path.circle(cent,2*rad)
    cut = matplotlib.path.Path(vertices=np.concatenate([outer.vertices, inner.vertices[::-1, ...]]),codes=np.concatenate([outer.codes, inner.codes]))
    patch = matplotlib.patches.PathPatch(cut, facecolor='#fff', edgecolor='black', )
    ax.add_patch(patch)

    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    
def set_pl_significances_params(fig):
    ax = fig.get_axes()[0]
    ax.set_title('Enrichments of compositional\nannotations in categories')
    ax.set_xticklabels(['      ','      ','      '])
    ax.set_xlabel('region')
    ax.set_yticklabels(['    ','    ','    '])
    ax.set_ylabel('type')
    ax.text(-0.20, 0.14+0/3, '\u2588', transform=ax.transAxes,color=colors_celltype[0], fontsize='x-large')
    ax.text(-0.20, 0.14+1/3, '\u2588', transform=ax.transAxes,color=colors_celltype[1], fontsize='x-large')
    ax.text(-0.20, 0.14+2/3, '\u2588', transform=ax.transAxes,color=colors_celltype[2], fontsize='x-large')
    ax.text(0.10+0/3, -0.20, '\u2588', transform=ax.transAxes,color=colors_region['0'], fontsize='x-large')
    ax.text(0.10+1/3, -0.20, '\u2588', transform=ax.transAxes,color=colors_region['1'], fontsize='x-large')
    ax.text(0.10+2/3, -0.20, '\u2588', transform=ax.transAxes,color=colors_region['2'], fontsize='x-large')
