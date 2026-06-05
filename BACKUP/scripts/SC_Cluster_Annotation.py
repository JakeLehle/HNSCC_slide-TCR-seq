#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 25 09:16:17 2025

@author: jlehle
"""

#%% Read in the SRAscraper config.yaml to get the working dir
# This section is the only thing that need to be automated to bring it into the snakmake pipeline
import os
import yaml
# This is the one line you would have to change vvvvvvv
os.chdir('/work/sdz852/WORKING/SC/fastq/Breast_cancer')

with open('config_SABCS.yaml', 'r') as config:
    config_yaml = yaml.safe_load(config)

#%% Go load the meatdata of all the files downloaded
working_dir = config_yaml['output_dir']
figures_dir = os.path.join(working_dir, 'figures')
os.makedirs(figures_dir, exist_ok=True)
os.chdir(working_dir)
os.getcwd()

import pickle
with open(os.path.join(working_dir, 'metadata/dictionary_file_filtered.pkl'), 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')
     
#%% Convert dict to dataframe
import pandas as pd

x_data = pd.DataFrame()
for key in gse_dict.keys():
    gse_dict[key]['series_id'] = key
    x_data = pd.concat([x_data, gse_dict[key]])
x_data.reset_index(inplace=True)
x_data = x_data.drop(columns=['index'])
x_data['run_alias'] = x_data['run_alias'].str.split('_').str[0]


#%% Here I will define the funtion to read in the regular single cell experiment results coming from cell ranger
from tqdm import tqdm
import scanpy as sc
import anndata as ad

def create_adata_normal_sc():
    failed_samples = []
    successful_samples = []
    
    # Step 1: Find and initialize with FIRST SUCCESSFUL sample
    adata = None
    first_sample_idx = None
    
    for i in range(len(x_data['series_id'])):
        data_file = os.path.join(working_dir, 'fastq', x_data['series_id'][i], 
                                 x_data['run_accession'][i], 
                                 x_data['run_accession'][i]+'_S1_L001_/outs/filtered_feature_bc_matrix')
        try:
            adata = sc.read_10x_mtx(data_file)
            
            # Add metadata
            for metadata in x_data.columns:
                adata.obs[metadata] = [x_data[metadata][i]]*adata.n_obs
            
            # Create unique obs names
            result = list(map("-".join, zip(adata.obs_names.to_list(), 
                                           adata.obs.run_accession.to_list())))
            adata.obs_names = result
            adata.var['gene_symbol'] = adata.var.index
            
            first_sample_idx = i
            successful_samples.append(i)
            tqdm.write(f"  Initialized with sample {i}: {x_data['run_accession'][i]}")
            break
            
        except Exception as e:
            tqdm.write(f"   Failed to read sample {i}: {x_data['run_accession'][i]}")
            tqdm.write(f"   Error: {str(e)}")
            failed_samples.append({
                'index': i,
                'run_accession': x_data['run_accession'][i],
                'series_id': x_data['series_id'][i],
                'error': str(e)
            })
    
    if adata is None:
        raise RuntimeError("No valid samples found to initialize AnnData object!")
    
    # Step 2: Loop through REMAINING samples
    for i in tqdm(range(first_sample_idx + 1, len(x_data['series_id'])), desc='Reading anndata'):
        data_file = os.path.join(working_dir, 'fastq', x_data['series_id'][i], 
                                 x_data['run_accession'][i], 
                                 x_data['run_accession'][i]+'_S1_L001_/outs/filtered_feature_bc_matrix')
        
        try:
            adata_tmp = sc.read_10x_mtx(data_file)
            
            # Add metadata
            for metadata in x_data.columns:
                adata_tmp.obs[metadata] = [x_data[metadata][i]]*adata_tmp.n_obs
            
            # Create unique obs names
            result = list(map("-".join, zip(adata_tmp.obs_names.to_list(), 
                                           adata_tmp.obs.run_accession.to_list())))
            adata_tmp.obs_names = result
            adata_tmp.var['gene_symbol'] = adata_tmp.var.index
            
            # Concatenate
            adata = ad.concat([adata, adata_tmp], join='outer', merge='same')
            
            successful_samples.append(i)
            
        except Exception as e:
            tqdm.write(f"   Failed to read sample {i}: {x_data['run_accession'][i]}")
            tqdm.write(f"   Error: {str(e)}")
            failed_samples.append({
                'index': i,
                'run_accession': x_data['run_accession'][i],
                'series_id': x_data['series_id'][i],
                'error': str(e)
            })
    
    # Filter gse_dict and save (same as before)
    x_data_filtered = x_data.iloc[successful_samples].reset_index(drop=True)
    
    successful_series_ids = set(x_data_filtered['series_id'].unique())
    gse_dict_filtered = {}
    
    for series_id in successful_series_ids:
        if series_id in gse_dict:
            successful_runs = set(x_data_filtered[x_data_filtered['series_id'] == series_id]['run_accession'])
            original_df = gse_dict[series_id]
            filtered_df = original_df[original_df['run_accession'].isin(successful_runs)].copy()
            
            if len(filtered_df) > 0:
                gse_dict_filtered[series_id] = filtered_df
    
    # Save filtered dictionary
    with open(os.path.join(working_dir, 'metadata/dictionary_file_filtered_post_adata.pkl'), 'wb') as pkl_file:
        pickle.dump(gse_dict_filtered, pkl_file)
    print(f"Filtered dictionary saved with {len(gse_dict_filtered)} series IDs")
    
    # Save metadata
    x_data_filtered.to_pickle(os.path.join(working_dir, 'metadata/sample_df_file_post_adata.pkl'))
    x_data_filtered.to_csv(os.path.join(working_dir, 'metadata/sample_df_file_post_adata.tsv'), 
                          sep='\t', index=False)
    
    # Save failed samples
    if failed_samples:
        failed_df = pd.DataFrame(failed_samples)
        failed_df.to_csv(os.path.join(working_dir, 'metadata/failed_samples.tsv'), 
                        sep='\t', index=False)
        print(f"\n{len(failed_samples)} samples failed. Details saved to metadata/failed_samples.tsv")
    
    print(f"\nSuccessfully loaded {len(successful_samples)} out of {len(x_data['series_id'])} samples")
    
    return adata, x_data_filtered, gse_dict_filtered
#%% Create the adata object and filtered dictionary
adata, x_data, gse_dict = create_adata_normal_sc()
#%% Save your work

for target in list(adata.obs.columns):
    adata.obs[target] = [str(element) for element in adata.obs[target]]

adata.obs.columns = adata.obs.columns.str.replace('/', '_')

adata.write("adata.h5ad")

#%% Let's write a function to do all the QC
import numpy as np
from scipy.stats import median_abs_deviation
import matplotlib.pyplot as plt
import seaborn as sb

def is_outlier(normal_anndata_object, metric: str, nmads: int):
    M = normal_anndata_object.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def QC_on_adata_normal(normal_anndata_object):
    normal_anndata_object.var["mt"] = normal_anndata_object.var_names.str.startswith(('mt-', 'MT-', 'Mt-'))
    # ribosomal genes
    normal_anndata_object.var["ribo"] = normal_anndata_object.var_names.str.startswith(('Rps', 'Rpl', 'RPS', 'RPL'))
    # hemoglobin genes.
    normal_anndata_object.var["hb"] = normal_anndata_object.var_names.str.contains(('^HB[^(P)]'))
    #We can now calculate the respective QC metrics with scanpy.
    sc.pp.calculate_qc_metrics(
        normal_anndata_object, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
    print(f"Generating figures with sample metrics saved in: {working_dir+'/figures'}")
    sc.pl.scatter(normal_anndata_object, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save='Non-norm_scatter_before_QC')
    plt.rcParams["figure.figsize"] = (20, 6)
    sc.pl.violin(normal_anndata_object, "pct_counts_mt", groupby='source_name', rotation=90, save='Non-norm_violin_pct_mito_before_QC')
    sc.pl.violin(normal_anndata_object, 'total_counts', groupby='source_name', rotation=90, save='Non-norm_violin_total_counts_before_QC')
    plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]
    p1 = sb.displot(normal_anndata_object.obs["total_counts"], bins=100, kde=False)
    p1.fig.savefig('Non-norm_hist_plot_total_counts_before_QC.pdf')
    os.rename(working_dir + '/Non-norm_hist_plot_total_counts_before_QC.pdf', working_dir + '/figures/Non-norm_hist_plot_total_counts_before_QC.pdf')
    # We now apply this function to the log1p_total_counts, log1p_n_genes_by_counts 
    # and pct_counts_in_top_20_genes QC covariates each with a threshold of 5 MADs.
    normal_anndata_object.obs["outlier"] = (
        is_outlier(normal_anndata_object, "log1p_total_counts", 5)
        | is_outlier(normal_anndata_object, "log1p_n_genes_by_counts", 5)
        | is_outlier(normal_anndata_object, "pct_counts_in_top_20_genes", 5)
    )
    normal_anndata_object.obs.outlier.value_counts()

    # pct_counts_Mt is filtered with 3 MADs. Additionally, cells with a percentage
    # of mitochondrial counts exceeding 20 % are filtered out.

    normal_anndata_object.obs["mt_outlier"] = is_outlier(normal_anndata_object, "pct_counts_mt", 3) | (
        normal_anndata_object.obs["pct_counts_mt"] > 20
    )
    normal_anndata_object.obs.mt_outlier.value_counts()

    # We now filter our AnnData object based on these two additional columns.
    print(f"Total number of cells: {normal_anndata_object.n_obs}")
    normal_anndata_object = normal_anndata_object[(~normal_anndata_object.obs.outlier) & (~normal_anndata_object.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {normal_anndata_object.n_obs}")

    # Final round of filtering
    sc.pp.filter_cells(normal_anndata_object, min_genes=200)
    sc.pp.filter_genes(normal_anndata_object, min_cells=20)
    print(f"Number of genes after cell filter: {normal_anndata_object.n_vars}")
    print(f"Final number of cells: {normal_anndata_object.n_obs}")


    sc.pl.scatter(normal_anndata_object, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save='Non-norm_scatter_after_QC')
    plt.rcParams["figure.figsize"] = (20, 6)
    sc.pl.violin(normal_anndata_object, 'pct_counts_mt', groupby='source_name', rotation=90, save='Non-norm_violin_pct_mito_after_QC')
    sc.pl.violin(normal_anndata_object, 'total_counts', groupby='source_name', rotation=90, save='Non-norm_violin_total_counts_after_QC')
    plt.rcParams["figure.figsize"] = plt.rcParamsDefault["figure.figsize"]
    p2 = sb.displot(normal_anndata_object.obs["total_counts"], bins=100, kde=False)
    p2.fig.savefig('Non-norm_hist_plot_total_counts_after_QC.pdf')
    os.rename(working_dir + '/Non-norm_hist_plot_total_counts_after_QC.pdf', working_dir + '/figures/Non-norm_hist_plot_total_counts_after_QC.pdf')
    
    return normal_anndata_object

#%% Run the QC on the adata object and normalize
adata_tmp = QC_on_adata_normal(adata)

#%% Automated cell type analysis with popv
import popv
# popV needs adata.X raw information no normalization
# use the QC_on_adata_normal() function to remove low quality cells as a pre-processing step 
# You will start with the adata_tmp object for this then

# Select a pre-trained model
huggingface_repo = "popV/tabula_sapiens_All_Cells"
# The query batch key is what will be used by bbknn for batch correction
query_batch_key = "run_accession"
algorithms = None
#%% Perform annotation useing a premade model
import numba
hmo = popv.hub.HubModel.pull_from_huggingface_hub(huggingface_repo, cache_dir="tmp/tabula_sapiens")
#%%
adata_tmp_an = hmo.annotate_data(
    adata_tmp,
    query_batch_key=query_batch_key,
    prediction_mode="inference",  # "fast" does not integrate reference and query.
    gene_symbols='feature_name' # "Uncomment if using gene symbols."
)

adata_tmp_an.write("adata_popv_an.h5ad")

#%% merge the popv annotations with the adata_tmp object

adata_tmp_an = adata_tmp_an[adata_tmp_an.obs["_dataset"] == "query"]
adata_tmp_an.obs = adata_tmp_an.obs.drop(columns=list(adata_tmp.obs.columns))
merged_df = pd.merge(adata_tmp.obs, adata_tmp_an.obs, left_index=True, right_index=True, how='inner')
adata_tmp = adata_tmp[adata_tmp.obs_names.isin(merged_df.index)]
merged_df = merged_df.loc[adata_tmp.obs_names]
adata_tmp.obs = merged_df
assert all(adata_tmp.obs_names == merged_df.index), "Indices do not match!"
#%% Pre-process samples to prep them for knn clustering
adata_pp = adata_tmp.copy()
sc.pp.normalize_total(adata_pp, target_sum=1e6)
sc.pp.log1p(adata_pp)

#%% Calculate the neighborhood plot and cluster
from os import cpu_count

NCPUS = cpu_count()
sc.pp.neighbors(adata_pp, n_pcs=NCPUS)
sc.tl.umap(adata_pp)
sc.tl.leiden(adata_pp, key_added='clusters', resolution=1, random_state=42)


def sort_by_substrings(list_of_strings, substring):
    with_substring = []
    without_substring = []
    for string in list_of_strings:
        if substring in string:
            with_substring.append(string)
        else:
            without_substring.append(string)

    return with_substring + without_substring

sorted_sources = sort_by_substrings(list(adata_pp.obs.source_name.unique()), 'norm')

sorted_sources
adata_pp.obs['source_name_sorted'] = pd.Categorical(
    values=adata_pp.obs.source_name, categories=sorted_sources, ordered=True
)
cmap = plt.get_cmap('turbo')

# Generate evenly spaced values equal to the value of groups in the dataset
values = np.linspace(0, 1, len(adata_pp.obs.run_accession.unique()))

# Get RGB values for each value in the colormap
color_list = [cmap(value) for value in values]

ncol = 3
nrow = 1
figsize = 4
wspace = 2
# Adapt figure size based on number of rows and columns and added space between them
# (e.g. wspace between columns)
fig, axs = plt.subplots(
    nrow, ncol, figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
)
plt.subplots_adjust(wspace=wspace)
sc.pl.umap(adata_pp, color=['source_name_sorted'], cmap=cmap, size=5, ax=axs[0], show=False)
sc.pl.umap(adata_pp, color=['clusters'], cmap=cmap, size=5, ax=axs[1], show=False)
sc.pl.umap(adata_pp, color=['run_accession'], palette=color_list, size=5, ax=axs[2], save='UMAP_before_batch_correction')

sc.external.pp.bbknn(adata_pp, batch_key='run_accession', n_pcs=NCPUS)  # running bbknn 1.3.6
sc.tl.umap(adata_pp)

fig, axs = plt.subplots(
    nrow, ncol, figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
)
plt.subplots_adjust(wspace=wspace)
sc.pl.umap(adata_pp, color=['source_name_sorted'], cmap=cmap, ax=axs[0], show=False)
sc.pl.umap(adata_pp, color=['clusters'], cmap=cmap, ax=axs[1],  show=False)
sc.pl.umap(adata_pp, color=['run_accession'], palette=color_list, ax=axs[2], save='UMAP_after_batch_correction')

#%% display the automated cell types
import matplotlib.patheffects as pe
from adjustText import adjust_text
import matplotlib.pyplot as plt

def gen_mpl_labels(
    adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None, color_by_group=False
):
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}

    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_umap"], axis=0)

    # Fill the text colors dictionary
    text_colors = {group: None for group in adata.obs[groupby].cat.categories}

    if color_by_group and groupby + "_colors" in adata.uns:
        for i, group in enumerate(adata.obs[groupby].cat.categories):
            if group in exclude:
                continue
            text_colors[group] = adata.uns[groupby + "_colors"][i]

    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, color=text_colors[k], **text_kwargs) for k, (x, y) in medians.items()
        ]
    else:
        texts = [ax.text(x=x, y=y, s=k, color=text_colors[k], **text_kwargs) for k, (x, y) in medians.items()]

    adjust_text(texts, **adjust_kwargs)


        
def nonoverlapping_UMAP(adata_obj, group_name):
    cmap = plt.get_cmap('turbo')
    # Generate evenly spaced values equal to the value of groups in the dataset
    adata_obj.obs[group_name] = adata_obj.obs[group_name].cat.remove_unused_categories()
    value_cat = pd.Categorical(adata_obj.obs[group_name])
    values = np.linspace(0, 1, len(value_cat.categories))
    # Get RGB values for each value in the colormap
    palette = [cmap(value) for value in values]
    
    # Combined path effects - white border first, then black border
    combined_effects = [
        pe.withStroke(linewidth=6, foreground="white"),  # Thick white border (outer)
        pe.withStroke(linewidth=1, foreground="black"),  # Thin black border (inner)
        pe.Normal()  # Original text color
    ]
    
    with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": 150, "figure.frameon": False}):
        ax = sc.pl.umap(adata_obj, color=group_name, show=False, legend_loc=None, frameon=False, size=5, palette=palette)
        gen_mpl_labels(
            adata_obj,
            group_name,
            exclude=("None",),
            ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
            text_kwargs=dict(
                fontsize=26, 
                path_effects=combined_effects  # Use the combined effects here
            ),
            color_by_group=True
        )
        fig = ax.get_figure()
        fig.tight_layout()
        plt.show()

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import stats

def assign_cluster_based_annotation(adata_obj):
    # Make a copy of the relevant columns
    df = adata_obj.obs[['clusters', 'popv_prediction', 'popv_prediction_score']].copy()
    
    # Ensure scores are numeric
    df['popv_prediction_score'] = pd.to_numeric(df['popv_prediction_score'], errors='coerce')
    
    # Min-Max normalization within each cluster
    def min_max_normalize(x):
        x_min = x.min()
        x_range = x.max() - x_min
        # Handle case where all scores are identical (avoid division by zero)
        if x_range == 0:
            return np.zeros(len(x))
        return (x - x_min) / x_range
    
    df['normalized_score'] = df.groupby('clusters')['popv_prediction_score'].transform(min_max_normalize)
    
    # LINEAR weighting (simple normalization to sum to 1 within each cluster)
    def linear_weights(x):
        total = x.sum()
        # Handle case where all scores are zero
        if total == 0:
            return np.ones(len(x)) / len(x)  # Equal weights
        return x / total
    
    df['weight'] = df.groupby('clusters')['normalized_score'].transform(linear_weights)
    
    # Calculate weighted score
    df['weighted_score'] = df['popv_prediction_score'] * df['weight']
    
    # Aggregate scores by cluster and predicted cell type
    cluster_type_scores = df.groupby(['clusters', 'popv_prediction'])['weighted_score'].sum().reset_index()
    
    # For each cluster, get the cell type with highest aggregated weighted score
    dominant_types = cluster_type_scores.loc[
        cluster_type_scores.groupby('clusters')['weighted_score'].idxmax()
    ].set_index('clusters')['popv_prediction']
    
    # Map the dominant type back to all cells
    adata_obj.obs['final_annotation'] = adata_obj.obs['clusters'].map(dominant_types)
    
    return adata_obj
#%%
# Apply the function to your data
adata_pp = assign_cluster_based_annotation(adata_pp)

nonoverlapping_UMAP(adata_pp, 'final_annotation')

adata_pp

#%% Save your work

#for target in list(adata_pp.obs.columns):
#    adata_pp.obs[target] = [str(element) for element in adata_pp.obs[target]]

adata_pp.write("adata_pp.h5ad")

#%% Map the final_annotation onto the adata_tmp object

import sys

adata_tmp = adata_tmp[adata_tmp.obs_names.isin(adata_pp.obs_names)]

sc.pp.normalize_total(adata_tmp, target_sum=1e6)

try:
    assert all(adata_tmp.obs_names == adata_pp.obs_names), "Indices do not match!"
except AssertionError as e:
    print(f"Assertion failed: {e}")
    # Exit the script with a non-zero status code indicating an error
    sys.exit(1)  
    

adata_tmp.obs['final_annotation'] = adata_pp.obs['final_annotation'].loc[adata_tmp.obs_names]

adata_tmp.write("adata_tmp.h5ad")

#%% Plotting section #%%

#%% plot for UMAP of popv annotation
cmap = plt.get_cmap('turbo')
value_cat = pd.Categorical(adata_pp.obs['popv_prediction'])
values = np.linspace(0, 1, len(value_cat.categories))
# Get RGB values for each value in the colormap
palette = [cmap(value) for value in values]

from matplotlib.patches import Patch
with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": 150, "figure.frameon": False}):
    fig = sc.pl.umap(adata_pp, color=['popv_prediction'], palette=palette, frameon=False, size=5, return_fig=True)
    
    plt.tight_layout()  # rect=[left, bottom, right, top] (0.85 leaves 15% space on right)
    plt.savefig(os.path.join(figures_dir, 'UMAP_popv_annotation.pdf'), bbox_inches='tight')
    plt.show()


#%% plot for UMAP of popv prediction score
cmap = plt.get_cmap('magma')
value_cat = pd.Categorical(adata_pp.obs['popv_prediction_score'])
values = np.linspace(0, 1, len(value_cat.categories))
# Get RGB values for each value in the colormap
palette = [cmap(value) for value in values]

from matplotlib.patches import Patch
with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": 150, "figure.frameon": False}):
    fig = sc.pl.umap(adata_pp, color=['popv_prediction_score'], palette=palette, frameon=False, size=5, legend_loc=None, return_fig=True)
    
    # Manually create legend with squares
    legend_elements = [
        Patch(facecolor='#000001', edgecolor='k', label='1'),
        Patch(facecolor='#230752', edgecolor='k', label='2'),
        Patch(facecolor='#6d1078', edgecolor='k', label='3'),
        Patch(facecolor='#ac2469', edgecolor='k', label='4'),
        Patch(facecolor='#d34737', edgecolor='k', label='5'),
        Patch(facecolor='#fd7c17', edgecolor='k', label='6'),
        Patch(facecolor='#ffbf14', edgecolor='k', label='7'),
        Patch(facecolor='#fafea0', edgecolor='k', label='8')
    ]
    
    ax = fig.axes[0]
    
    # Place legend to the right of the plot
    legend = ax.legend(
        handles=legend_elements,
        fontsize=26,
        frameon=False,
        bbox_to_anchor=(1, 1),  # Coordinates for legend position (1.05 = right outside)
        loc='upper left',           # Anchor point for the legend
        borderaxespad=0.            # Padding between legend and axes
    )
    
    # Adjust the figure to make room for the legend
    plt.tight_layout(rect=[0, 0, 1.2, 1])  # rect=[left, bottom, right, top] (0.85 leaves 15% space on right)
    plt.savefig(os.path.join(figures_dir, 'UMAP_popv_prediction_score.pdf'), bbox_inches='tight')
    plt.show()
#%% Pre process which clusters have the most cells per each final annotation to make stacked bar plot
import pandas as pd
from collections import defaultdict

def get_dominant_clusters(adata_pp):
    """
    For each unique cell type in 'final_annotation', finds the cluster 
    (from 'clusters') that contains the largest number of cells of that type.
    
    Returns:
        List of cluster IDs with the dominant population for each cell type
    """
    # Create a DataFrame with just the annotations and clusters
    df = pd.DataFrame({
        'cell_type': adata_pp.obs['final_annotation'],
        'cluster': adata_pp.obs['clusters']
    })
    
    # Group by cell type and cluster, then count occurrences
    cluster_counts = df.groupby(['cell_type', 'cluster']).size().reset_index(name='count')
    
    # For each cell type, find the cluster with maximum count
    dominant_clusters = (
        cluster_counts
        .sort_values('count', ascending=False)
        .drop_duplicates('cell_type')
        ['cluster']
        .tolist()
    )
    
    # Convert to strings and remove duplicates while preserving order
    seen = set()
    unique_dominant_clusters = [
        str(x) for x in dominant_clusters 
        if not (str(x) in seen or seen.add(str(x)))
    ]
    
    return unique_dominant_clusters

# Usage:
CLUSTERS = get_dominant_clusters(adata_pp)
print(f"Dominant clusters for each cell type: {CLUSTERS}")
#%% Stacked bar plots to validate final annotation cell type dominance excluded weighting scores (presentational)

import pandas as pd
from matplotlib.colors import to_hex
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap

MIN_PERCENT = 5  # Label threshold
LABEL_SPACING = 10  # Vertical spacing between labels
BAR_WIDTH = 0.6
PLOT_SPACING = 1.0  # Space between cluster bars

# Create figure
fig, ax = plt.subplots(figsize=(24, 8))  # Width adjusted for 12 clusters

# Get color mapping
all_celltypes = adata_pp.obs['popv_prediction'].astype('category').cat.categories
cmap = plt.get_cmap('turbo')
values = np.linspace(0, 1, len(all_celltypes))
color_map = dict(zip(all_celltypes, [to_hex(cmap(value)) for value in values]))

for cluster_idx, cluster_num in enumerate(CLUSTERS):
    # Filter data for current cluster
    sub_clust_cells = adata_pp.obs[adata_pp.obs['clusters'] == cluster_num]
    celltype_counts = sub_clust_cells['popv_prediction'].value_counts()
    percentages = (celltype_counts / celltype_counts.sum() * 100).round(2)
    
    # Process labels and colors
    nonzero_mask = percentages > 0
    filtered_labels = percentages.index[nonzero_mask].tolist()
    filtered_percentages = percentages[nonzero_mask].tolist()
    filtered_colors = [color_map[label] for label in filtered_labels]

    # Filter for significant segments
    significant_idx = [i for i,p in enumerate(filtered_percentages) if p >= MIN_PERCENT]
    labels = [filtered_labels[i] for i in significant_idx]
    percentages_sig = [filtered_percentages[i] for i in significant_idx]
    colors_sig = [filtered_colors[i] for i in significant_idx]

    # Sort largest to smallest
    sorted_idx = np.argsort(percentages_sig)[::-1]
    labels = [labels[i] for i in sorted_idx]
    percentages_sig = [percentages_sig[i] for i in sorted_idx]
    colors_sig = [colors_sig[i] for i in sorted_idx]

    # X-position for this cluster
    x_pos = cluster_idx * PLOT_SPACING
    
    # Draw stacked bar
    bottom = 0
    for percent, color in zip(filtered_percentages, filtered_colors):
        alpha = 0.4 if percent < MIN_PERCENT else 1.0
        ax.bar(x_pos, percent, width=BAR_WIDTH, bottom=bottom, 
               color=color, edgecolor='white', alpha=alpha)
        bottom += percent

    # Calculate segment centers
    segment_centers = np.cumsum(filtered_percentages) - np.array(filtered_percentages)/2

    # Trace types annotation at bottom center
    small_total = sum(p for p in filtered_percentages if p < MIN_PERCENT)
    if small_total > 0:
        ax.text(x_pos, 3, 
                f"Trace cells \ntotal: {small_total:.1f}%",
                va='bottom', ha='center', fontsize=18 , style='italic',
                bbox=dict(facecolor='white', alpha=0.9, edgecolor='lightgray'))

    # Place labels centered in the bar
    for idx, label in zip(significant_idx, labels):
        percent = filtered_percentages[idx]
        color = filtered_colors[idx]
        y_center = segment_centers[idx]
        
        # Only show label if there's enough space
        if percent >= MIN_PERCENT:  # Only label segments at least twice the minimum
            wrapped_label = '\n'.join(wrap(label, width=12))  # More aggressive wrapping
            ax.text(x_pos, y_center, 
                    f"{wrapped_label}\n({percent:.1f}%)",
                    va='center', ha='center', fontsize=20 if percent > 20 else 14,
                    color='black',
                    bbox=dict(facecolor="white", alpha=0.9,
                             edgecolor='lightgray', boxstyle='round,pad=0.2'))


    # Cluster label below bars
    ax.text(x_pos, -5, f"{cluster_num}", 
            ha='center', va='top', fontsize=36)

# Final formatting
ax.set_xlim(-0.5, len(CLUSTERS)*PLOT_SPACING - (PLOT_SPACING-BAR_WIDTH))
ax.set_ylim(0, 100)
ax.set_ylabel('Percentage (%)', fontsize=36)
ax.spines[['top', 'right']].set_visible(False)
ax.grid(axis='y', linestyle=':', alpha=0.3)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

plt.tight_layout()
plt.savefig(os.path.join(figures_dir, 'Stacked_Bar_Plot_Support_For_Final_Annotation.pdf'), 
            bbox_inches='tight', dpi=300)
plt.show()
#%% End file
import sys

sys.exit()
#%%
# #%%
# print("Loading processed AnnData object...")
adata_pp = sc.read_h5ad(filename=os.path.join(working_dir, 'adata_pp.h5ad'))
nonoverlapping_UMAP(adata_pp, 'final_annotation')
adata_pp.obs['final_annotation'].value_counts()
# adata_pp.obs['final_annotation'] = adata_pp.obs['final_annotation'].str.replace('ovarian surface epithelial cell', 'surface epithelial cell')
