import warnings
warnings.filterwarnings('ignore')

import dynamo as dyn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

test_adata_tmp = dyn.read_h5ad('label.h5ad')
test_adata_tmp.obs['label_time'] = 2 # this is the labeling time 
tkey = 'label_time'
cell_filter = (test_adata_tmp.X> 0).sum(1) > 800
gene_list = pd.read_csv('total_gene_list.txt', header=None)
test_adata_tmp_DE = test_adata_tmp[cell_filter.A1, gene_list]
###### new
dyn.pp.recipe_monocle(test_adata_tmp_DE, tkey='label_time', experiment_type='one-shot',reset_X=True)
dyn.tl.dynamics(test_adata_tmp_DE, model='auto',cores=1,NTR_vel = True)
dyn.tl.reduceDimension(test_adata_tmp_DE, reduction_method='umap')
dyn.tl.cell_velocities(test_adata_tmp_DE, calc_rnd_vel=True)
test_adata_tmp_DE.obs['time'] = test_adata_tmp_DE.obs.time.astype('category')
dyn.pl.umap(test_adata_tmp_DE, color='time', save_show_or_return='show', color_key_cmap = 'viridis')
dyn.pl.streamline_plot(test_adata_tmp_DE, color='time', color_key_cmap = 'viridis', 
                       basis='umap', show_legend='right', save_show_or_return='return')
dyn.pl.streamline_plot(test_adata_tmp_DE, color='time', color_key_cmap = 'viridis', basis='umap_rnd',
                       show_legend='right', save_show_or_return='return')
dyn.tl.cell_wise_confidence(test_adata_tmp_DE, ekey='M_t', vkey='velocity_N', method='correlation')
dyn.pl.umap(test_adata_tmp_DE, color='correlation_velocity_confidence',pointsize=0.8)
new_gene_list = pd.read_csv("new_gene_list.txt",header=None)
dyn.pl.phase_portraits(test_adata_tmp_DE, genes=list(new_gene_list), color='time',  discrete_continous_div_color_key=[None, None, None], 
                       discrete_continous_div_color_key_cmap=['viridis', None, None], ncols=6,  pointsize=5,
                       save_show_or_return='return')

