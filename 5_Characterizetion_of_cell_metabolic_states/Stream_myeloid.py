import stream as st
# st.__version__
import os 
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

import sys
firboblast_file=sys.argv[0]
work_dir=sys.argv[1]

st.set_figure_params(dpi=80,style='white',figsize=[5.4,3],
                     rc={'image.cmap': 'viridis'})
os.chdir(firboblast_file)
adata = sc.read_loom("./myeloid_integration_subset.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene')
st.add_cell_labels(adata,file_name='./myeloid_integration.celltype.subset.tsv')
st.set_workdir(adata,work_dir)
adata.X=adata.X.toarray()

loess_frac=0.05
n_genes=[1500,2000,2500,3000,3500,4000,4500,5000]
n_pc=[10,15,20,25,30,35,40,45,50]
n_neighbors=[50,60,70,80,90,100]
n_components=4
n_clusters=12

def get_trajectory(loess_frac,n_genes,n_pc,n_neighbors,n_components,n_clusters):
    st.select_variable_genes(adata,loess_frac=loess_frac,n_genes=n_genes)
    st.select_top_principal_components(adata,feature='var_genes',n_pc=n_pc,first_pc=True)
    st.dimension_reduction(adata,method='se',feature='var_genes',n_neighbors=n_neighbors, n_components=n_components,n_jobs=50)
    st.plot_dimension_reduction(adata,color=['label'],show_graph=False,show_text=False)
    st.plot_visualization_2D(adata,method='umap',n_neighbors=n_neighbors,color=['label'],use_precomputed=False)
    st.plot_visualization_2D(adata,method='umap',n_neighbors=n_neighbors,color=['label'],use_precomputed=False)
    st.seed_elastic_principal_graph(adata,n_clusters=n_clusters,use_vis=True)
    st.plot_dimension_reduction(adata,color=['label'],n_components=2,show_graph=True,show_text=False)
    st.plot_branches(adata,show_text=True)
    st.elastic_principal_graph(adata,epg_alpha=0.015,epg_mu=0.1,epg_lambda=0.01)
    st.plot_dimension_reduction(adata,color=['label'],n_components=2,show_graph=True,show_text=False)
    st.plot_branches(adata,show_text=True)
    st.add_cell_labels(adata,file_name='./myeloid_integration.celltype.subset.tsv')
    adata.uns['label_color']['Mono_CD16']='#92C274'
    adata.uns['label_color']['Mono_FCN1']='#1C7E76'
    adata.uns['label_color']['Mono_THBS1']='#824DAE'
    adata.uns['label_color']['Macro_SPP1']='#DB5DB8'
    adata.uns['label_color']['Macro_SLPI']='#CA131F'
    adata.uns['label_color']['Macro_C1QC']='#8ACDD4'
    adata.uns['label_color']['Macro_IL1B']='#1965B0'
    adata.uns['label_color']['Macro_IL32']='#FBB065'
    adata.uns['label_color']['Macro_CDC20']='#FED43B'
    adata.uns['label_color']['Macro_CLEC10A']='#B8A1CD'
    adata.uns['label_color']['Macro_CCL18']='#FD8586'
    st.plot_stream_sc(adata,root='S0',color=['label'],
                      dist_scale=0.5,show_graph=True,show_text=True)    
    st.plot_stream(adata,root='S0',color=['label'])
    plt.savefig(work_dir,'myeloid_stream_celltype_'+
                str(loess_frac)+'_'+str(n_genes)+'_'+str(n_pc)+'_'+str(n_neighbors)+
                '_'+str(n_components)+'_'+str(n_clusters)+'.pdf', format='pdf', bbox_inches='tight')
    st.write(adata,file_name=work_dir+myeloid_stream_metatype_'+
            str(loess_frac)+'_'+str(n_genes)+'_'+str(n_pc)+'_'+str(n_neighbors)+
            '_'+str(n_components)+'_'+str(n_clusters)+'.pkl')
    adata.obs.to_csv(work_dir+'myeloid_stream_metatype_'+
                str(loess_frac)+'_'+str(n_genes)+'_'+str(n_pc)+'_'+str(n_neighbors)+
                '_'+str(n_components)+'_'+str(n_clusters)+'.cells_transition.csv', index=True)
    
for gene in n_genes:
    for pc in n_pc:
        for neighbors in n_neighbors:
            file_name = work_dir+'myeloid_stream_metatype_' + \
                        str(loess_frac) + '_' + str(gene) + '_' + str(pc) + '_' + str(neighbors) + \
                        '_' + str(n_components) + '_' + str(n_clusters) + '.pkl'
            if not os.path.exists(file_name):  # check files
                get_trajectory(loess_frac, gene, pc, neighbors, n_components, n_clusters)
                
## use the parameters
loess_frac=0.05
n_genes=3500
n_pc=20
n_neighbors=80
n_components=4
n_clusters=12

st.add_cell_labels(adata,file_name='./myeloid_integration.celltype.subset.tsv')
adata.uns['label_color']['Mono_CD16']='#92C274'
adata.uns['label_color']['Mono_FCN1']='#1C7E76'
adata.uns['label_color']['Mono_THBS1']='#824DAE'
adata.uns['label_color']['Macro_SPP1']='#DB5DB8'
adata.uns['label_color']['Macro_SLPI']='#CA131F'
adata.uns['label_color']['Macro_C1QC']='#8ACDD4'
adata.uns['label_color']['Macro_IL1B']='#1965B0'
adata.uns['label_color']['Macro_IL32']='#FBB065'
adata.uns['label_color']['Macro_CDC20']='#FED43B'

st.plot_stream(adata,root='S0',color=['label'])
plt.savefig(work_dir+'/myeloid_stream_celltype_'+
            str(loess_frac)+'_'+str(n_genes)+'_'+str(n_pc)+'_'+str(n_neighbors)+
            '_'+str(n_components)+'_'+str(n_clusters)+'.pdf', format='pdf', bbox_inches='tight')

st.add_cell_labels(adata,file_name='./myeloid_integration.metatype.subset.tsv')
adata.uns['label_color']['PUFA_LTM']='#1C7E76'
adata.uns['label_color']['PUFA_EM']='#92C274'
adata.uns['label_color']['OXP']='Thistle'
adata.uns['label_color']['FAO']='#824DAE'
adata.uns['label_color']['SLM']='#8ACDD4'
adata.uns['label_color']['ARG']='#DB5DB8'
adata.uns['label_color']['FAS']='#608FBF'
adata.uns['label_color']['GST']='#CA131F'
adata.uns['label_color']['PUFA']='#FBB065'
adata.uns['label_color']['AA']='#FED43B'

st.plot_stream(adata,root='S1',color=['label'])
plt.savefig(work_dir+'myeloid_stream_metatype_'+
            str(loess_frac)+'_'+str(n_genes)+'_'+str(n_pc)+'_'+str(n_neighbors)+
            '_'+str(n_components)+'_'+str(n_clusters)+'.pdf', format='pdf', bbox_inches='tight')

st.plot_stream(adata,root='S1',color=['S1_pseudotime'],factor_min_win=1.5)
plt.savefig(work_dir+'stream_psedotime_'+'S1_pseudotime'+'.pdf', format='pdf', bbox_inches='tight')

new_features=pd.read_csv(firboblast_file+'myeloid_integration_sub.gsva.csv',
                        index_col=0)
adata.obs = pd.concat([adata.obs, new_features], axis=1)