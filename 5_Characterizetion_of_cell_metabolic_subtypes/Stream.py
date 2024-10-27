# Analyzes the trajectory of cell lineages using STREAM.
#' @param cell_lineage The cell lineage to analysis.
#' @param file Path to the cell lineages object file.
#' @param work_dir Path to the output directory.
#' @author Ke Tang
#
import stream as st
# st.__version__
import os 
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

import sys

cell_lineage=sys.argv[0]
file=sys.argv[1]
work_dir=sys.argv[2]

st.set_figure_params(dpi=80,style='white',figsize=[5.4,3],rc={'image.cmap': 'viridis'})
os.chdir(file)
adata = sc.read_loom(str('./'+cell_lineage+"_integration.loom"), sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene')
st.add_cell_labels(adata,file_name=str('./'+cell_lineage+'_integration.celltype.tsv'))
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
    st.dimension_reduction(adata,method='se',feature='var_genes',n_neighbors=n_neighbors, n_components=n_components,n_jobs=10)
    st.plot_dimension_reduction(adata,color=['label'],show_graph=False,show_text=False)
    st.plot_visualization_2D(adata,method='umap',n_neighbors=n_neighbors,color=['label'],use_precomputed=False)
    st.plot_visualization_2D(adata,method='umap',n_neighbors=n_neighbors,color=['label'],use_precomputed=False)
    st.seed_elastic_principal_graph(adata,n_clusters=n_clusters,use_vis=True)
    st.plot_dimension_reduction(adata,color=['label'],n_components=2,show_graph=True,show_text=False)
    st.plot_branches(adata,show_text=True)
    st.elastic_principal_graph(adata,epg_alpha=0.015,epg_mu=0.1,epg_lambda=0.01)
    st.plot_dimension_reduction(adata,color=['label'],n_components=2,show_graph=True,show_text=False)
    st.plot_branches(adata,show_text=True)
    st.add_cell_labels(adata,file_name=str('./'+cell_lineage+'_integration.celltype.tsv'))
    st.plot_stream_sc(adata,root='S0',color=['label'],
                      dist_scale=0.5,show_graph=True,show_text=True)    
    st.plot_stream(adata,root='S0',color=['label'])
    plt.savefig(work_dir+str('./'+cell_lineage+'_stream_celltype_')+
                str(loess_frac)+'_'+str(n_genes)+'_'+str(n_pc)+'_'+str(n_neighbors)+
                '_'+str(n_components)+'_'+str(n_clusters)+'.pdf', format='pdf', bbox_inches='tight')
    st.write(adata,file_name=work_dir+str('./'+cell_lineage+'_stream_metatype_')+
            str(loess_frac)+'_'+str(n_genes)+'_'+str(n_pc)+'_'+str(n_neighbors)+
            '_'+str(n_components)+'_'+str(n_clusters)+'.pkl')
    adata.obs.to_csv(work_dir+str('./'+cell_lineage+'_stream_metatype_')+
                str(loess_frac)+'_'+str(n_genes)+'_'+str(n_pc)+'_'+str(n_neighbors)+
                '_'+str(n_components)+'_'+str(n_clusters)+'.cells_transition.csv', index=True)
    
for gene in n_genes:
    for pc in n_pc:
        for neighbors in n_neighbors:
            file_name = work_dir+str('./'+cell_lineage+'_stream_metatype_') + \
                        str(loess_frac) + '_' + str(gene) + '_' + str(pc) + '_' + str(neighbors) + \
                        '_' + str(n_components) + '_' + str(n_clusters) + '.pkl'
            if not os.path.exists(file_name):  # check files
                get_trajectory(loess_frac, gene, pc, neighbors, n_components, n_clusters)