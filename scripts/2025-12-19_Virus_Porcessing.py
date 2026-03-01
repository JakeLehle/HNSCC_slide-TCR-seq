#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 19 21:14:55 2025

@author: sdz852
"""
import os
import yaml
# Load config and setup paths
config_path = '/work/sdz852/WORKING/SC/fastq/Breast_cancer/config_SABCS.yaml'
with open(config_path, 'r') as config:
    config_yaml = yaml.safe_load(config)
working_dir = config_yaml['output_dir']
os.chdir(working_dir)
#%% Go load the meatdata of all the files downloaded
figures_dir = os.path.join(working_dir, 'figures')
os.makedirs(figures_dir, exist_ok=True)
os.chdir(working_dir)
os.getcwd()
import pickle
with open(os.path.join(working_dir, 'metadata/dictionary_file_filtered_post_adata.pkl'), 'rb') as pkl_file:
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

# Clean column names to be HDF5-compatible
print("Cleaning metadata column names for HDF5 compatibility...")
rename_map = {}
for col in x_data.columns:
    # Replace forward slashes, backslashes, and clean up spaces
    new_col = str(col).replace('/', '_').replace('\\', '_').replace(' ', '_')
    # Remove other potentially problematic characters
    new_col = new_col.replace('(', '').replace(')', '').replace('[', '').replace(']', '')
    if new_col != col:
        rename_map[col] = new_col
        print(f"  Renamed: '{col}' -> '{new_col}'")

if rename_map:
    x_data.rename(columns=rename_map, inplace=True)
    print(f"Cleaned {len(rename_map)} column names")
else:
    print("No problematic column names found")
#%% Prep for reading in all the kraken2 files
# read in the hierarchy.txt file with all he viruses used form the kraken2 database. This is output in each folder and is the same for all of them so I just need to read in one.
with open(os.path.join(working_dir, 'fastq', x_data['series_id'][0], x_data['run_accession'][0], x_data['run_accession'][0]+'_S1_L001_/outs/kraken2_filtered_feature_bc_matrix/hierarchy.txt'), 'rt') as f:
    lines = f.readlines()
import re
column_4 = []
column_5 = []
for line in lines:
    # Assuming your columns are separated by spaces
    #row = re.split(r'\t+', line.rstrip('\n'))
    pattern = re.compile('[^\t]+')
    pattern_2 = re.compile(r'.*?([a-zA-Z].*)')
    row = pattern.findall(line.rstrip('\n'))
    # Replace 0 with the column index you want to extract (starting from 0)
    column_4.append(row[4])
    column_5.append(pattern_2.findall(row[5])[0])
     
v_hierarchy = {key: value for key, value in zip(column_4, column_5)}
#%% Make sure that all files have all viruses if they don't then you wont be able to concatenate them easily. 
counter = 0
import gzip
import csv

# read the tsv
for key in x_data['series_id'].unique():
    for accession in x_data['run_accession'][x_data['series_id'] == key]:
        try:
            genes_tsv_path = os.path.join(working_dir, 'fastq', key, accession, 
                                        accession+"_S1_L001_/outs/kraken2_filtered_feature_bc_matrix/genes.tsv.gz")
            matrix_mtx_path = os.path.join(working_dir, 'fastq', key, accession, 
                                         accession+"_S1_L001_/outs/kraken2_filtered_feature_bc_matrix/matrix.mtx.gz")
            
            if not os.path.exists(genes_tsv_path):
                raise FileNotFoundError(f"genes.tsv.gz not found for {accession}")
            if not os.path.exists(matrix_mtx_path):
                raise FileNotFoundError(f"matrix.mtx.gz not found for {accession}")

            accession_list = []
            with gzip.open(genes_tsv_path, 'rt') as f:
                tsv_reader = csv.reader(f, delimiter="\t")
                for row in tsv_reader:
                    accession_list.append(row[1])
                accession_list_final = [x for x in accession_list if x.startswith(tuple(v_hierarchy.values()))]
                
                for tax_id, virus in v_hierarchy.items():
                    if virus in accession_list_final:
                        print(f'{virus} reads found in {accession}')
                    else:
                        with gzip.open(genes_tsv_path, 'at') as f:
                            tsv_writer = csv.writer(f, delimiter="\t", lineterminator="\n")
                            tsv_writer.writerow([tax_id, virus])       
            
            with gzip.open(genes_tsv_path, 'rt') as f:
                feature_tsv = []
                tsv_reader = csv.reader(f, delimiter="\t")
                for row in tsv_reader:
                    feature_tsv.append(row[1])
                feature_tsv_len = len(feature_tsv)
                print(f'Finished writing {accession} : gene.tsv.gz')
            
            with gzip.open(matrix_mtx_path) as file:
                lines = file.readlines()
                tmp = str(feature_tsv_len) + ' ' + str(lines[3 - 1]).split(' ', 2)[1] + ' ' + str(len(lines) - 3) + '\n'
                lines[3 - 1] = bytes(tmp, 'utf-8')
            
            with gzip.open(matrix_mtx_path, "wb") as file:
                for line in lines:
                    file.write(line)
                print(f'Finished writing {accession} : matrix.mtx.gz')
            
            counter = counter + 1
            print(f'Counter: {counter}')
            
        except FileNotFoundError as e:
            print(f"ERROR: Skipping {accession} - {str(e)}")
            continue
        except Exception as e:
            print(f"ERROR: Unexpected error processing {accession}: {str(e)}")
            continue 
#%% Make the function that will create the adata_v object
import gzip
import csv
from tqdm import tqdm
import anndata as ad
import scanpy as sc
import os

def create_adata_viral_sc():
    adata_v = None
    
    for i in tqdm(range(len(x_data['series_id'])), desc='Reading viral anndata'):
        try:
            # Build paths and check files
            data_file = os.path.join(working_dir, 'fastq', x_data['series_id'][i], 
                                   x_data['run_accession'][i], 
                                   x_data['run_accession'][i]+'_S1_L001_/outs/kraken2_filtered_feature_bc_matrix')
            
            # Check if directory exists
            if not os.path.exists(data_file):
                print(f"ERROR: Directory not found - {data_file}")
                continue
                
            file_list = ['barcodes.tsv.gz', 'genes.tsv.gz', 'matrix.mtx.gz']
            
            # Check all required files exist
            missing_files = [f for f in file_list if not os.path.exists(os.path.join(data_file, f))]
            if missing_files:
                print(f"ERROR: Missing files in {data_file}: {', '.join(missing_files)}")
                continue
            
            # Process files
            for file in file_list:
                try:
                    with gzip.open(os.path.join(data_file, file), 'rt') as f_in, \
                         open(os.path.join(data_file, file[:-3]), 'wt') as f_out:
                        f_out.writelines(f_in)
                except Exception as e:
                    print(f"ERROR: Failed to decompress {file} in {data_file}: {str(e)}")
                    # Clean up partially extracted files
                    if os.path.exists(os.path.join(data_file, file[:-3])):
                        os.remove(os.path.join(data_file, file[:-3]))
                    continue
            
            # Read the data
            try:
                adata_v_tmp = sc.read_10x_mtx(data_file)
                
                # Add metadata
                for metadata in x_data.columns:
                    adata_v_tmp.obs[metadata] = [x_data[metadata][i]]*adata_v_tmp.n_obs
                
                # Process observation names
                result = list(map("-".join, zip(adata_v_tmp.obs_names.to_list(), 
                                             adata_v_tmp.obs.run_accession.to_list())))
                adata_v_tmp.obs_names = result    
                adata_v_tmp.var['gene_symbol'] = adata_v_tmp.var.index
                
                # Concatenate with main object
                if adata_v is None:
                    adata_v = adata_v_tmp
                else:
                    adata_v = ad.concat([adata_v, adata_v_tmp], join='outer', merge='same')
                
            except Exception as e:
                print(f"ERROR: Failed to process 10x data in {data_file}: {str(e)}")
                continue
            
            # Clean up extracted files
            for file in file_list:
                extracted_file = os.path.join(data_file, file[:-3])
                if os.path.exists(extracted_file):
                    os.remove(extracted_file)
                    
        except Exception as e:
            print(f"ERROR: Unexpected error processing sample {x_data['run_accession'][i]}: {str(e)}")
            continue
    
    if adata_v is None:
        print("ERROR: No valid data was processed - returning empty AnnData object")
        return ad.AnnData()
    
    return adata_v
#%% Create the adata_v object
adata_v = create_adata_viral_sc()

#%% Save your work
for target in list(adata_v.obs.columns):
    adata_v.obs[target] = [str(element) for element in adata_v.obs[target]]

adata_v.write("adata_v.h5ad")

#%% load back in with all datasets
import scanpy as sc
adata = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata.h5ad')
)

adata_v = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_v.h5ad')
)

adata_pp = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_pp_CD.h5ad')
)
#%% Find all human viruses using the human kraken database
# This will need to be an input provided by the user on the ClusterCatcher pipeline. This will be a herarchy file from kraken2 of human sepecific viruses.  
with open('/work/sdz852/WORKING/kraken2/human_viral/inspect.txt', 'rt') as f:
    lines = f.readlines()

import re
v_hierarchy_human = {}
v_hierarchy_human_species_list = []

for line in lines:
    # Skip empty lines
    if not line.strip():
        continue
        
    # Split on whitespace (but keep the last columns together)
    parts = line.strip().split()
    
    # We need at least 6 columns to process
    if len(parts) < 6:
        continue
        
    # Extract columns (adjust indices based on actual structure)
    # The pattern appears to be: percentage, read_count, tax_count, rank_code, tax_id, name
    percentage = parts[0]
    read_count = parts[1]
    tax_count = parts[2]
    rank_code = parts[3]  # This is the S/D/K/etc. column
    tax_id = parts[4]
    # The name is everything after tax_id (columns 5+)
    name = ' '.join(parts[5:])
    
    # Store in dictionary with tax_id as key
    v_hierarchy_human[tax_id] = {
        'rank_code': rank_code,
        'name': name,
        'read_count': read_count,
        'tax_count': tax_count,
        'percentage': percentage
    }
    
    # Add to species list if it's a species (starts with S)
    if rank_code.startswith('S'):
        v_hierarchy_human_species_list.append({
            'tax_id': tax_id,
            'name': name,
            'read_count': read_count
        })

# Create the name lists you originally wanted
v_hierarchy_human_list = tuple(entry['name'] for entry in v_hierarchy_human.values())
v_hierarchy_human_species_names = tuple(entry['name'] for entry in v_hierarchy_human_species_list)

# I'll also make a list of all the gene names from adata_pp which will come in handy later
adata_pp_gene_names = list(adata_pp.var.index)

#%% Strip all viruses not found in humans away from the adata_v object
human_virus_names = tuple(name.strip() for name in v_hierarchy_human_species_names if name.strip())
# Create a mask for human viruses in adata_v
human_virus_mask = adata_v.var_names.isin(human_virus_names)
# Drop all the non-human viruses
adata_v = adata_v[:, human_virus_mask].copy()
#%% Revert adata_pp log normalization to get raw counts
import anndata as ad
import numpy as np

# Revert log1p transformation: exp(x) - 1
adata_pp.X = np.expm1(adata_pp.X)

# Note: This gives you normalized counts (e.g., CPM or counts per 10k)
# If you need truly raw counts, you'd need to reverse the normalize_total step too
# by multiplying by the normalization factor, but for joining with viral data
# and re-normalizing together, this intermediate state is actually ideal

print(f"adata_pp shape: {adata_pp.shape}")
print(f"adata_pp data range: {adata_pp.X.min():.2f} to {adata_pp.X.max():.2f}")

#%% Prepare viral data (keep as raw counts - DO NOT normalize yet)
# Filter cells to match adata_pp
adata_v_filtered = adata_v[adata_v.obs.index.isin(list(adata_pp.obs_names))].copy()

# Filter to only cells with viral reads (optional - you may want to keep all cells)
sc.pp.filter_cells(adata_v_filtered, min_counts=1)

print(f"adata_v_filtered shape: {adata_v_filtered.shape}")
print(f"adata_v_filtered data range: {adata_v_filtered.X.min():.2f} to {adata_v_filtered.X.max():.2f}")

#%% Filter adata_pp by cell barcodes present in adata_v_filtered
adata_filtered = adata_pp[adata_pp.obs.index.isin(list(adata_v_filtered.obs_names))].copy()

# Clear any varm data that might cause concatenation issues
adata_filtered.varm = {}
adata_v_filtered.varm = {}

print(f"adata_filtered shape: {adata_filtered.shape}")
print(f"Cell barcode overlap: {len(set(adata_filtered.obs_names) & set(adata_v_filtered.obs_names))}")

#%% Join the unnormalized data
adata_joined = ad.concat(
    [adata_filtered, adata_v_filtered],
    join='outer',
    axis=1,
    merge='same'
)

# Reorder to match adata_pp
common_cells = list(set(adata_pp.obs_names) & set(adata_joined.obs_names))
adata_joined = adata_joined[common_cells, :].copy()

# Transfer metadata from adata_pp
adata_joined.obs = adata_pp.obs.loc[common_cells].copy()

assert all(adata_joined.obs_names == adata_pp.obs.loc[common_cells].index), (
    f"Order mismatch after reordering!"
)

print(f"adata_joined shape after merging: {adata_joined.shape}")
print(f"adata_joined data range: {adata_joined.X.min():.2f} to {adata_joined.X.max():.2f}")

#%% Now normalize the joined data together
print("Normalizing joined dataset...")

# Normalize total counts per cell
sc.pp.normalize_total(adata_joined, target_sum=1e4)

# Log transform
sc.pp.log1p(adata_joined)

print(f"After normalization - data range: {adata_joined.X.min():.2f} to {adata_joined.X.max():.2f}")
#%% Calculate the neighborhood plot and cluster
from os import cpu_count
import matplotlib.pyplot as plt

NCPUS = cpu_count()
sc.pp.neighbors(adata_joined, n_pcs=NCPUS)
sc.tl.umap(adata_joined)
sc.tl.leiden(adata_joined, key_added='clusters', resolution=1, random_state=42)

def sort_by_substrings(list_of_strings, substring):
    with_substring = []
    without_substring = []
    for string in list_of_strings:
        if substring in string:
            with_substring.append(string)
        else:
            without_substring.append(string)
    return with_substring + without_substring

sorted_sources = sort_by_substrings(list(adata_joined.obs.source_name.unique()), 'norm')

adata_joined.obs['source_name_sorted'] = pd.Categorical(
    values=adata_joined.obs.source_name, categories=sorted_sources, ordered=True
)
cmap = plt.get_cmap('turbo')

# Generate evenly spaced values equal to the value of groups in the dataset
values = np.linspace(0, 1, len(adata_joined.obs.run_accession.unique()))
color_list = [cmap(value) for value in values]

ncol = 3
nrow = 1
figsize = 4
wspace = 2

# BEFORE batch correction
fig, axs = plt.subplots(
    nrow, ncol, figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
)
plt.subplots_adjust(wspace=wspace)
sc.pl.umap(adata_joined, color=['source_name_sorted'], cmap=cmap, size=5, ax=axs[0], show=False)
sc.pl.umap(adata_joined, color=['clusters'], cmap=cmap, size=5, ax=axs[1], show=False)
sc.pl.umap(adata_joined, color=['run_accession'], palette=color_list, size=5, ax=axs[2], show=False)
fig.savefig(os.path.join(figures_dir, 'UMAP_before_batch_correction.png'), dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {os.path.join(figures_dir, 'UMAP_before_batch_correction.png')}")

# Apply batch correction
sc.external.pp.bbknn(adata_joined, batch_key='run_accession')
sc.tl.umap(adata_joined)

# AFTER batch correction
fig, axs = plt.subplots(
    nrow, ncol, figsize=(ncol * figsize + (ncol - 1) * wspace * figsize, nrow * figsize)
)
plt.subplots_adjust(wspace=wspace)
sc.pl.umap(adata_joined, color=['source_name_sorted'], cmap=cmap, ax=axs[0], show=False)
sc.pl.umap(adata_joined, color=['clusters'], cmap=cmap, ax=axs[1], show=False)
sc.pl.umap(adata_joined, color=['run_accession'], palette=color_list, ax=axs[2], show=False)
fig.savefig(os.path.join(figures_dir, 'UMAP_after_batch_correction.png'), dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {os.path.join(figures_dir, 'UMAP_after_batch_correction.png')}")

#%% Rank the genes
sc.tl.rank_genes_groups(adata_joined, groupby='final_annotation', key_added='rank_genes', method='wilcoxon')

for target in list(adata_joined.obs.columns):
    adata_joined.obs[target] = [str(element) for element in adata_joined.obs[target]]

if 'gene_ids' in adata_joined.var.columns:
    adata_joined.var['gene_ids'] = [str(item) for item in adata_joined.var['gene_ids']]

adata_joined.write("adata_v_pp.h5ad")
print(f"Saved: adata_v_pp.h5ad")

#%% Calculate the differential expression of the viruses across multiple cell types
from tqdm import tqdm
import numpy as np
results = adata_joined.uns['rank_genes']
remove_list = adata_pp_gene_names
out = np.array([[0,0,0,0,0]])

for group in tqdm(results['names'].dtype.names, desc="Processing DEGs"):
    out = np.vstack((out, np.vstack((results['names'][group],
                                     results['scores'][group],
                                     results['pvals_adj'][group],
                                     results['logfoldchanges'][group],
                                     np.array([group] * len(results['names'][group])).astype('object'))).T))

markers = pd.DataFrame(out[1:], columns = ['virus', 'scores', 'pval_adj', 'lfc', 'cluster'])
markers = markers[(markers.scores > 0.0001)]
mask = markers.virus.str.contains('|'.join(remove_list)) 
markers = markers[~mask]
markers['pval_adj'] = markers['pval_adj'].map('{:e}'.format)
adata_joined.uns['markers'] = markers

#%% Filter for virus markers and check if any found
import sys
virus_markers = markers[markers['virus'].isin(human_virus_names)]

if len(virus_markers) == 0:
    print("\n" + "="*80)
    print("WARNING: No human viruses detected in dataset")
    print("="*80)
    # Save empty results
    summary_file = os.path.join(working_dir, 'virus_detection_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("No human viruses detected in this dataset\n")
        f.write(f"Total markers analyzed: {len(markers)}\n")
        f.write(f"Human virus database size: {len(human_virus_names)}\n")
    print(f"Saved summary: {summary_file}")
    sys.exit(0)

# Continue with virus analysis
subset = virus_markers.iloc[:, [0,1]]
subset = subset.groupby('virus')['scores'].sum()

print(f"\nTop 10 viruses detected:")
for i, virus in enumerate(list(subset.sort_values(ascending=False).index)[:10], 1):
    print(f"  {i}. {virus}: {subset[virus]:.4f}")

# Save virus detection summary
summary_file = os.path.join(working_dir, 'virus_detection_summary.txt')
with open(summary_file, 'w') as f:
    f.write(f"Virus Detection Summary\n")
    f.write(f"="*80 + "\n\n")
    f.write(f"Total viruses detected: {len(subset)}\n")
    f.write(f"Total cells analyzed: {adata_joined.n_obs}\n\n")
    f.write("Top 10 viruses by aggregate score:\n")
    for i, virus in enumerate(list(subset.sort_values(ascending=False).index)[:10], 1):
        f.write(f"  {i}. {virus}: {subset[virus]:.4f}\n")
print(f"Saved summary: {summary_file}")

#%% Plotting section to display virus
adata_joined.var.index = adata_joined.var['gene_symbol']

# Matrix plot - Top 10 viruses
plt.rcParams['figure.figsize'] = (8, 8)
sc.settings.verbosity = 3
sc.set_figure_params(scanpy=True, fontsize=16)

fig = sc.pl.matrixplot(
    adata_joined, 
    list(subset.sort_values(ascending=False).index)[:10], 
    'final_annotation', 
    dendrogram=True, 
    var_group_rotation=30, 
    cmap='plasma', 
    log=True,
    return_fig=True,
    show=False
)
matrix_plot_path = os.path.join(figures_dir, 'Matrix_Plot_Top_10_Virus.png')
fig.savefig(matrix_plot_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {matrix_plot_path}")

# Top virus UMAP
Top_virus_gene = list(subset.sort_values(ascending=False).index)[0]
Top_virus_counts = adata_joined[:, Top_virus_gene].X.toarray().flatten()
Top_virus_series = pd.Series(Top_virus_counts, index=adata_joined.obs.index, name='Top_virus_counts')

adata_pp.obs[Top_virus_gene] = adata_pp.obs.index.map(Top_virus_series).fillna(0)
print(f"\nTop virus: {Top_virus_gene}")
print(f"Added {Top_virus_gene} counts to adata_pp.obs")
print(f"Non-zero cells: {(adata_pp.obs[Top_virus_gene] > 0).sum()}")
print(f"Stats for positive cells:\n{adata_pp.obs[Top_virus_gene][adata_pp.obs[Top_virus_gene] > 0].describe()}")

sc.set_figure_params(scanpy=True, fontsize=25)
fig = sc.pl.umap(
    adata_joined, 
    color=Top_virus_gene, 
    use_raw=False, 
    size=5, 
    color_map='plasma', 
    frameon=False,
    return_fig=True,
    show=False
)
umap_virus_path = os.path.join(figures_dir, f'UMAP_{Top_virus_gene.replace(" ", "_")}.png')
fig.savefig(umap_virus_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {umap_virus_path}")

plt.rcdefaults()

#%% Save final processed object
adata_pp.write(os.path.join(working_dir, "adata_pp_with_virus.h5ad"))
print(f"\nSaved final object: adata_pp_with_virus.h5ad")

#%% End of the script
print("\n" + "="*80)
print("Pipeline completed successfully!")
print("="*80)
import sys 
sys.exit(0)