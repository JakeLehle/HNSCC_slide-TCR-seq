#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 09:48:01 2025

@author: sdz852
"""
#!/usr/bin/env python
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
    
#%% load back in with all datasets
import scanpy as sc
adata = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata.h5ad')
)

adata_tmp = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_tmp.h5ad')
)

adata_pp = sc.read_h5ad(
    filename=os.path.join(working_dir, 'adata_pp_with_virus.h5ad')
)

#%% Do some QC to check how this looks on our single cell data
# display the automated cell types
import matplotlib.patheffects as pe
from adjustText import adjust_text
import matplotlib.pyplot as plt
import numpy as np
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
    cmap = plt.get_cmap('tab20b')
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
    
    with plt.rc_context({"figure.figsize": (10, 10), "figure.dpi": 150, "figure.frameon": False}):
        ax = sc.pl.umap(adata_obj, color=group_name, show=False, legend_loc=None, title='', frameon=False, size=5, palette=palette)
        gen_mpl_labels(
            adata_obj,
            group_name,
            exclude=("None",),
            ax=ax,
            adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
            text_kwargs=dict(
                fontsize=36, 
                path_effects=combined_effects  # Use the combined effects here
            ),
            color_by_group=True
        )
        fig = ax.get_figure()
        fig.tight_layout()
        plt.savefig(os.path.join(figures_dir, 'Final_Cell_Type_Annotation_UMAP.pdf'))
        plt.show()
        plt.close()

#%% COSMIC mutation signature section 
import pandas as pd
import numpy as np
from collections import defaultdict

def reverse_complement(sequence):
    """Generate reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in sequence[::-1]])
    
def get_mutation_type(ref, alt, context):
    """
    Determine the mutation type in the format required by SigProfiler
    Format: [5' base][REF>ALT][3' base]
    """
    if len(ref) != 1 or len(alt) != 1:
        return None
    
    # For pyrimidine contexts (C or T)
    if ref in ['C', 'T']:
        return f"{context[0]}[{ref}>{alt}]{context[2]}"
    # For purine contexts (A or G) - we should have converted these earlier
    else:
        return None

def process_mutations(input_file, output_file):
    # First read just the header line
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                # Remove the '#' and strip whitespace
                header = line.lstrip('#').strip().split('\t')
                break
    
    # Now read the data with the proper column names
    df = pd.read_csv(input_file, sep='\t', comment='#', names=header)
    
    # Filter out indels (keep only single nucleotide variants)
    df = df[(df['REF'].str.len() == 1) & (df['ALT_expected'].str.len() == 1)]
    # Filter out rows with 'N' in REF_TRI or ALT_TRI
    print(f"Total rows after indel filtering: {len(df)}")
    df = df[~df['REF_TRI'].str.contains('N') & ~df['ALT_TRI'].str.contains('N')]
    print(f"Rows after filtering out 'N' in trinucleotide contexts: {len(df)}")
    
    # Convert G or A mutations to their reverse complements
    for idx, row in df.iterrows():
        ref = row['REF']
        alt = row['ALT_expected']
        
        # Check if mutation needs to be converted (purine to pyrimidine)
        if ref in ['G', 'A']:
            # Convert reference and alternate alleles
            df.at[idx, 'REF'] = reverse_complement(ref)
            df.at[idx, 'ALT_expected'] = reverse_complement(alt)
            
            # Convert the trinucleotide context
            df.at[idx, 'REF_TRI'] = reverse_complement(row['REF_TRI'])
            df.at[idx, 'ALT_TRI'] = reverse_complement(row['ALT_TRI'])
    
    # Process each mutation to get the 96 mutation types
    mutation_counts = defaultdict(lambda: defaultdict(int))
    
    for _, row in df.iterrows():
        cb = row['CB']  # Cell barcode
        ref = row['REF']
        alt = row['ALT_expected']
        context = row['REF_TRI']  # Reference trinucleotide context
        
        # Skip if not a valid substitution
        if ref not in ['C', 'T']:
            continue
            
        # Skip if not one of the 6 possible mutation types
        if f"{ref}>{alt}" not in ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']:
            continue
            
        # Get the mutation type in SigProfiler format
        mutation_type = get_mutation_type(ref, alt, context)
        if mutation_type:
            mutation_counts[cb][mutation_type] += 1
    
    # Generate all possible 96 mutation types
    all_mutation_types = [
        f"{five}[{ref}>{alt}]{three}"
        for ref in ['C', 'T']
        for alt in (['A', 'G', 'T'] if ref == 'C' else ['A', 'C', 'G'])
        for five in ['A', 'C', 'G', 'T']
        for three in ['A', 'C', 'G', 'T']
    ]
    
    # Get all unique cell barcodes
    cell_barcodes = list(mutation_counts.keys())
    
    # Create a dictionary to hold our data
    data_dict = {mut_type: [] for mut_type in all_mutation_types}
    
    # Initialize counts to 0 for all mutation types
    for mut_type in all_mutation_types:
        data_dict[mut_type] = [mutation_counts[cb].get(mut_type, 0) for cb in cell_barcodes]
    
    # Create the DataFrame in one operation
    result_df = pd.DataFrame(data_dict, index=cell_barcodes).T  # Transpose to get mutations as rows
    
    # Ensure all mutation types are included (in case some weren't observed)
    result_df = result_df.reindex(all_mutation_types, fill_value=0)
    
    # Convert all values to integers
    result_df = result_df.astype(int)
    
    # Save to output file
    result_df.to_csv(output_file, sep='\t')
    
    return result_df
#%% Running the processing_mutations fucntion to make the sigprofiler_mutation_counts.txt
if __name__ == "__main__":
    input_file = "/work/sdz852/WORKING/SC/fastq/Breast_cancer/results_SABCS/all_samples.single_cell_genotype.filtered.tsv"
    output_file = "sigprofiler_mutation_counts.txt"
    
    print(f"Processing mutations from {input_file}...")
    result = process_mutations(input_file, output_file)
    print(f"Results saved to {output_file}")
    print(f"Final matrix shape: {result.shape}")
#%% Generate unique mutations table (removing cell barcode information)
import pandas as pd
import os

# Read in the original mutation file
input_file = "/work/sdz852/WORKING/SC/fastq/Breast_cancer/results_SABCS/all_samples.single_cell_genotype.filtered.tsv"

print(f"Reading mutation file: {input_file}")

# Read the file with proper header handling
with open(input_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            # Remove the '#' and strip whitespace
            header = line.lstrip('#').strip().split('\t')
            break

# Read the data with the proper column names
df = pd.read_csv(input_file, sep='\t', comment='#', names=header)

print(f"Original data shape: {df.shape}")
print(f"Columns: {df.columns.tolist()}")

df_no_cb = df.drop(columns=['CB', 'Num_cells_expected', 'Cell_type_observed', 'Base_observed', 'Num_reads'])

# Drop duplicate rows (now that CB is removed, duplicates represent the same mutation across cells)
df_unique = df_no_cb.drop_duplicates()

print(f"Shape after removing duplicates: {df_unique.shape}")
print(f"Number of duplicate rows removed: {df_no_cb.shape[0] - df_unique.shape[0]}")

# Save the unique mutations table
output_file = os.path.join(working_dir, "all_cell_types.single_cell_genotype.filtered.tsv")
df_unique.to_csv(output_file, sep='\t', index=False, header=True)

print(f"\nUnique mutations table saved to: {output_file}")
print(f"Final unique mutations: {df_unique.shape[0]}")

#%% Load the cell annotations file
cell_annotations_path = os.path.join(working_dir, 'cell_annotations.txt')
cell_annotations = pd.read_csv(cell_annotations_path, sep='\t')

# Extract the complete list of cell barcodes
all_expected_barcodes = set(cell_annotations['cell_barcodes'])

print(f"Total expected cells: {len(all_expected_barcodes)}")
#%% Confirm the sigprofiler_mutation_counts file was made and move forward with processing mutation matrix
result = pd.read_csv(os.path.join(working_dir, "sigprofiler_mutation_counts.txt"), sep= '\t')

mutations_per_cell = result.sum(axis=0)

mutations_per_cell_df = mutations_per_cell.to_frame(name='Total_Mutations')

# Make sure the cell barcodes match
mutations_per_cell_df = mutations_per_cell_df.reindex(adata_pp.obs.index.tolist())
mutations_per_cell_df['Total_Mutations'] = mutations_per_cell_df['Total_Mutations'].fillna(0)

# Add to AnnData obs
adata_pp.obs['total_mutations'] = mutations_per_cell_df['Total_Mutations']

#%% Load the combined callable sites file
callable_sites_path = os.path.join(working_dir, 'CombinedCallableSites', 'complete_callable_sites.tsv')
callable_sites_df = pd.read_csv(callable_sites_path, sep='\t')

# Get the cell barcodes from callable sites
callable_barcodes = set(callable_sites_df['CB'])

print(f"Cells in callable sites: {len(callable_barcodes)}")

# Find missing cells
len(list(set(callable_barcodes).intersection(set(mutations_per_cell[1:].index))))
len(list(set(mutations_per_cell[1:].index) - set(callable_barcodes)))

missing_cells = all_expected_barcodes - callable_barcodes
print(f"Cells missing from callable sites: {len(missing_cells)}")
missing_cells
# Check if these cells are in the mutation data

missing_intersection = list(set(missing_cells).intersection(set(mutations_per_cell[1:].index)))
print(f"Missing cells also absent from mutations: {len(missing_intersection)}")
len(set(missing_intersection).intersection(list(set(mutations_per_cell[1:].index) - set(callable_barcodes))))

#%% Okay I found two ways to get out barcodes where we have mutations recorded 
# but the cells don't meet the cutoff that is used for the count base function. 
# So this will be an issue and I need to go and modify the snp call function to meet the same strengency that a cell has to have at least 5x coverage at the 
# Before I increase the stringency I want to see how the data looks for the SNP calls as is. In the fututre I'm going to drop any cells where mutations were found but those mutations sites weren't > 5X in coverage

cells_to_drop = list(set(mutations_per_cell[1:].index) - set(callable_barcodes))

adata_pp.obs.loc[adata_pp.obs.index.isin(cells_to_drop), 'total_mutations'] = 0

# Identify columns to zero out
cols_to_zero = result.columns[1:].intersection(cells_to_drop)
result.loc[:, cols_to_zero] = 0
result.index = result['Unnamed: 0']
result = result.drop(columns="Unnamed: 0")

result.index.names = ['']

result.to_csv(os.path.join(working_dir, 'SNP_matrix_for_SigProfiler.txt'), sep='\t', index=True, header=True)
#%% Check to see if you had any changes at this point you should have 0 if eveything is working
from pandas.testing import assert_frame_equal
# Read back in the file you just made 
result_2 = pd.read_csv(os.path.join(working_dir, "SNP_matrix_for_SigProfiler.txt"), sep= '\t')
result_2.index = result_2['Unnamed: 0']
result_2 = result_2.drop(columns="Unnamed: 0")
result_2.index.names = ['']
assert_frame_equal(result, result_2)

# If the assert function threw an error stop and try to figure out what is going on
# as long as there was no differences you can use the results object if needed you can use the results_2 object to keep pushing the pipeline
callable_sites_df.index = callable_sites_df['CB']
callable_sites_df = callable_sites_df.reindex(adata_pp.obs.index.tolist())
callable_sites_df['CB'] = list(callable_sites_df.index)
callable_sites_df['SitesPerCell'] = callable_sites_df['SitesPerCell'].fillna(0)

#%% Plot some UMAPs **** CHECK THIS AND UPDATE ****
# Go ahead and add the total mutations to the adata for both the raw and then the normalized mutation counts 
adata_pp.obs['normalized_total_mutations'] = adata_pp.obs['total_mutations'] / callable_sites_df['SitesPerCell']
adata_pp.obs['normalized_total_mutations'] = adata_pp.obs['normalized_total_mutations'].fillna(0)
adata_pp.obs['normalized_total_mutations'] = np.nan_to_num(adata_pp.obs['normalized_total_mutations'], posinf=0.0, neginf=0.0)

sc.set_figure_params(scanpy=True, fontsize=25)
fig = sc.pl.umap(
    adata_pp, 
    color='total_mutations', 
    use_raw=False, 
    size=5, 
    color_map='plasma', 
    frameon=False,
    return_fig=True,
    show=False
)
total_mutation_path = os.path.join(figures_dir, 'UMAP_Total_Mutations.png')
fig.savefig(total_mutation_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {total_mutation_path}")


sc.set_figure_params(scanpy=True, fontsize=25)
fig = sc.pl.umap(
    adata_pp, 
    color='normalized_total_mutations', 
    use_raw=False, 
    size=5, 
    color_map='plasma', 
    frameon=False,
    return_fig=True,
    show=False
)
normalized_total_mutation_path = os.path.join(figures_dir, 'UMAP_Normalized_Total_Mutations.png')
fig.savefig(normalized_total_mutation_path, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {normalized_total_mutation_path}")
#%% Mutation Count Distribution Analysis for NMF Preprocessing
"""
This script analyzes your sparse mutation matrix to help identify cells
with sufficient mutation counts for effective NMF signature extraction.

Key outputs:
1. Distribution statistics of mutations per cell
2. Breakdown by mutation count thresholds
3. Sparsity analysis at different thresholds
4. Visualizations to guide threshold selection
5. Filtered matrix export for NMF analysis

Run interactively in your Python IDE (Spyder, PyCharm, VSCode, etc.)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set better plotting defaults
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

#%% Configuration
# This file should be set as an output of the SComatic pipeline and read in at the top of the script as an imput that is required
MATRIX_FILE = os.path.join(working_dir, "ALL_Cell_SNP_matrix_for_SigProfiler.txt")
# This looks like it is the same as the working dir right now
OUTPUT_DIR = working_dir
#%% Load Data
print("="*80)
print("LOADING MUTATION MATRIX")
print("="*80)

# Load the mutation matrix
data = pd.read_csv(MATRIX_FILE, sep='\t', index_col=0)
print(f"\nOriginal matrix shape: {data.shape}")
print(f"Mutation contexts (rows): {data.shape[0]}")
print(f"Cells/samples (columns): {data.shape[1]}")

# Calculate Mutations Per Cell
print("\n" + "="*80)
print("CALCULATING MUTATIONS PER CELL")
print("="*80)

mutations_per_cell = data.sum(axis=0)
#%%
"""
Semi-Supervised Signature Refitting for Single-Cell Data

Instead of de novo signature extraction, this approach:
1. Uses KNOWN COSMIC signatures relevant to HNSCC
2. Solves for per-cell weights (exposures) using NNLS
3. Evaluates how well H × W approximates the original mutation matrix X

This version includes:
- Scree plot-based signature selection (elbow detection)
- Fixed SBS40 signature detection
- Proper evaluation using Eckart-Young theorem

Author: Jake Lehle
Date: 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
from scipy.optimize import nnls
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, r2_score
import warnings
warnings.filterwarnings('ignore')

# COSMIC signature colors
COSMIC_COLORS = {
    'C>A': '#1EBFF0',
    'C>G': '#050708',
    'C>T': '#E62725',
    'T>A': '#CBCACB',
    'T>C': '#A1CE63',
    'T>G': '#EDB6C2'
}


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def log(msg, level="INFO"):
    """Print timestamped log message"""
    timestamp = datetime.now().strftime('%H:%M:%S')
    print(f"[{timestamp}] {level}: {msg}", flush=True)


def save_figure(fig, figure_dir, filename, dpi=300):
    """Save figure to output directory"""
    filepath = Path(figure_dir) / filename
    fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
    log(f"Saved: {filepath}")
    return filepath


# =============================================================================
# STEP 1: EXTRACT HNSCC-RELEVANT COSMIC SIGNATURES
# =============================================================================

def extract_hnscc_signatures(cosmic_file, output_dir=None):
    """
    Extract HNSCC-relevant COSMIC signatures
    
    HNSCC signatures based on literature:
    SBS1, SBS2, SBS3, SBS4, SBS5, SBS7a, SBS7b, SBS7d, SBS13, 
    SBS16, SBS17a, SBS17b, SBS18, SBS33, SBS40
    
    NOTE: SBS40 may appear as 'SBS40a', 'SBS40b', 'SBS40c' in COSMIC v3.4
    
    Parameters:
    -----------
    cosmic_file : str
        Path to COSMIC_v3.4_SBS_GRCh38.txt
    output_dir : str, optional
        Directory to save extracted signatures
    
    Returns:
    --------
    pd.DataFrame : HNSCC signature matrix (96 contexts × signatures)
    """
    log("="*80)
    log("EXTRACTING HNSCC-RELEVANT COSMIC SIGNATURES")
    log("="*80)
    
    # HNSCC-relevant signatures
    hnscc_signatures = [
        "SBS1", "SBS2", "SBS3", "SBS4", "SBS5", 
        "SBS6", "SBS7a", "SBS7b", "SBS7c", "SBS7d",
        "SBS8", "SBS9", "SBS10a", "SBS10b", "SBS10c", 
        "SBS10d", "SBS11", "SBS12", "SBS13", "SBS14", 
        "SBS15", "SBS16", "SBS17a", "SBS17b", "SBS18", 
        "SBS19", "SBS20", "SBS21", "SBS22a", "SBS22b", 
        "SBS23", "SBS24", "SBS25", "SBS26", 
        "SBS28", "SBS29", "SBS30", "SBS31", "SBS32", 
        "SBS33", "SBS34", "SBS35", "SBS36", "SBS37", 
        "SBS38", "SBS39", "SBS40a", "SBS40b", "SBS40c", 
        "SBS41", "SBS42", "SBS44", 
        "SBS84", "SBS85", "SBS86", "SBS87", "SBS88", 
        "SBS89", "SBS90", "SBS91", "SBS92", "SBS93", 
        "SBS94", "SBS96", "SBS97", "SBS98", 
        "SBS99", "SBS100", "SBS101", "SBS102", "SBS103", 
        "SBS104", "SBS105", "SBS106", "SBS107", "SBS108", 
        "SBS109", "SBS110"]

    log(f"\nLoading COSMIC signatures from: {cosmic_file}")
    cosmic_all = pd.read_csv(cosmic_file, sep='\t', index_col=0)
    log(f"Total COSMIC signatures available: {cosmic_all.shape[1]}")
    
    # Check for SBS40 variants (SBS40a, SBS40b, SBS40c)
    available_hnscc = []
    missing_hnscc = []
    
    for sig in hnscc_signatures:
        if sig in cosmic_all.columns:
            available_hnscc.append(sig)
        elif sig == 'SBS40':
            # Check for SBS40 variants
            variants = [col for col in cosmic_all.columns if col.startswith('SBS40')]
            if variants:
                log(f"\nℹ️  SBS40 variants found: {', '.join(variants)}", level="INFO")
                # Use SBS40a by default (most common)
                if 'SBS40a' in variants:
                    available_hnscc.append('SBS40a')
                    log(f"Using SBS40a as SBS40 representative")
                else:
                    available_hnscc.append(variants[0])
                    log(f"Using {variants[0]} as SBS40 representative")
            else:
                missing_hnscc.append(sig)
        else:
            missing_hnscc.append(sig)
    
    log(f"\nHNSCC signatures found: {len(available_hnscc)}/{len(hnscc_signatures)}")
    
    if available_hnscc:
        log(f"Available: {', '.join(available_hnscc)}")
    
    if missing_hnscc:
        log(f"⚠️  Missing: {', '.join(missing_hnscc)}", level="WARNING")
    
    # Extract HNSCC signatures
    hnscc_sigs = cosmic_all[available_hnscc].copy()
    
    log(f"\nExtracted signature matrix: {hnscc_sigs.shape[0]} contexts × {hnscc_sigs.shape[1]} signatures")
    
    # Verify normalization (should sum to 1)
    col_sums = hnscc_sigs.sum(axis=0)
    if not np.allclose(col_sums, 1.0, atol=1e-5):
        log("⚠️  Signatures not normalized, normalizing now...", level="WARNING")
        hnscc_sigs = hnscc_sigs / col_sums
    
    log("✓ Signatures are properly normalized (sum to 1.0)")
    
    # Save if output directory provided
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
        
        hnscc_file = output_path / "hnscc_cosmic_signatures.txt"
        hnscc_sigs.to_csv(hnscc_file, sep='\t', float_format='%.6f')
        log(f"\nSaved HNSCC signatures: {hnscc_file}")
        
        # Save list of signatures
        sig_list_file = output_path / "hnscc_signature_list.txt"
        with open(sig_list_file, 'w') as f:
            f.write("HNSCC-Relevant COSMIC Signatures\n")
            f.write("="*50 + "\n\n")
            for i, sig in enumerate(available_hnscc, 1):
                f.write(f"{i:2d}. {sig}\n")
        log(f"Saved signature list: {sig_list_file}")
    
    log("\n" + "="*80)
    log("HNSCC SIGNATURE EXTRACTION COMPLETE")
    log("="*80 + "\n")
    
    return hnscc_sigs


# =============================================================================
# STEP 2: SIGNATURE SELECTION VIA SCREE PLOT ELBOW DETECTION
# =============================================================================

def select_signatures_via_scree_plot(mutation_matrix, hnscc_sigs, core_signatures,
                                     candidate_pool, output_dir, max_signatures=15,
                                     verbose=True):
    """
    Select optimal number of signatures using scree plot elbow detection
    
    Strategy:
    1. Start with core signatures (e.g., SBS2, SBS3, SBS5)
    2. For each additional signature, calculate explained variance
    3. Find elbow point where adding more signatures gives diminishing returns
    4. Select best candidates based on individual explanatory power
    
    The key insight: More signatures ALWAYS reduce error, but we want to stop
    when the improvement becomes marginal (elbow point).
    
    Parameters:
    -----------
    mutation_matrix : pd.DataFrame
        Original mutation matrix X
    hnscc_sigs : pd.DataFrame
        All available HNSCC signatures
    core_signatures : list
        Core signatures to always include
    candidate_pool : list
        Pool of candidate signatures to consider
    output_dir : str
        Directory to save results
    max_signatures : int
        Maximum signatures to test (default: 15)
    verbose : bool
        Print detailed progress
    
    Returns:
    --------
    dict with selected signatures and analysis results
    """
    if verbose:
        log("="*80)
        log("SIGNATURE SELECTION VIA SCREE PLOT ELBOW DETECTION")
        log("="*80)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Verify signatures available
    available_sigs = set(hnscc_sigs.columns)
    core_available = [sig for sig in core_signatures if sig in available_sigs]
    candidates_available = [sig for sig in candidate_pool if sig in available_sigs and sig not in core_available]
    
    if verbose:
        log(f"\nSignature pool:")
        log(f"  Core signatures: {core_available}")
        log(f"  Candidate pool: {len(candidates_available)} signatures")
        log(f"  Max signatures to test: {max_signatures}")
    
    # -------------------------------------------------------------------------
    # Step 1: Rank candidates by individual explanatory power
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"STEP 1: Ranking candidates by individual explanatory power")
        log(f"{'='*70}")
    
    # For each candidate, fit it alone and see how much variance it explains
    candidate_scores = []
    
    X = mutation_matrix.values
    X_norm = np.linalg.norm(X, 'fro')
    
    for candidate_sig in candidates_available:
        # Fit just this signature
        sig_matrix = hnscc_sigs[[candidate_sig]]
        fitting = fit_signatures_nnls(mutation_matrix, sig_matrix, verbose=False)
        
        # Calculate explained variance
        residual = X - fitting['reconstruction'].values
        residual_norm = np.linalg.norm(residual, 'fro')
        explained_var = 1 - (residual_norm / X_norm)**2
        
        candidate_scores.append({
            'signature': candidate_sig,
            'explained_variance': explained_var,
            'residual_norm': residual_norm
        })
    
    # Sort by explained variance
    candidate_scores = sorted(candidate_scores, key=lambda x: x['explained_variance'], reverse=True)
    
    if verbose:
        log(f"\nCandidate ranking (by individual explained variance):")
        for i, score in enumerate(candidate_scores[:10], 1):
            log(f"  {i:2d}. {score['signature']:<10} explains {100*score['explained_variance']:.2f}% variance")
    
    # Create ordered candidate list
    candidates_ordered = [s['signature'] for s in candidate_scores]
    
    # -------------------------------------------------------------------------
    # Step 2: Build scree plot - test increasing signature counts
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"STEP 2: Building scree plot (testing 1 to {max_signatures} signatures)")
        log(f"{'='*70}")
    
    scree_data = []
    
    # Always start with core
    current_signatures = core_available.copy()
    
    # Add candidates one by one in order of explanatory power
    all_signatures_ordered = core_available + candidates_ordered
    
    for n_sigs in range(len(core_available), min(max_signatures + 1, len(all_signatures_ordered) + 1)):
        # Get signatures to test
        test_signatures = all_signatures_ordered[:n_sigs]
        sig_matrix = hnscc_sigs[test_signatures]
        
        # Fit and evaluate
        fitting = fit_signatures_nnls(mutation_matrix, sig_matrix, verbose=False)
        evaluation = evaluate_reconstruction(mutation_matrix, fitting['reconstruction'], verbose=False)
        
        # Calculate explained variance
        explained_var = 1 - evaluation['relative_frobenius_error']**2
        
        scree_data.append({
            'n_signatures': n_sigs,
            'signatures': test_signatures.copy(),
            'frobenius_error': evaluation['frobenius_error'],
            'relative_error': evaluation['relative_frobenius_error'],
            'explained_variance': explained_var,
            'optimality_ratio': evaluation['optimality_ratio']
        })
        
        if verbose and n_sigs <= 10:
            log(f"  {n_sigs:2d} signatures: Error={evaluation['frobenius_error']:>8.2f}, "
                f"RelErr={100*evaluation['relative_frobenius_error']:>5.2f}%, "
                f"ExpVar={100*explained_var:>5.2f}%")
    
    # -------------------------------------------------------------------------
    # Step 3: Find elbow point using second derivative
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"STEP 3: Detecting elbow point in scree plot")
        log(f"{'='*70}")
    
    errors = np.array([s['frobenius_error'] for s in scree_data])
    n_sigs_array = np.array([s['n_signatures'] for s in scree_data])
    
    # Calculate first derivative (rate of error decrease)
    dy = np.gradient(errors)
    
    # Calculate second derivative (rate of change of decrease)
    d2y = np.gradient(dy)
    
    # Find elbow: where second derivative is most positive (curvature changes most)
    # We want the point where error reduction starts slowing down significantly
    elbow_idx = np.argmax(d2y)
    elbow_n_sigs = scree_data[elbow_idx]['n_signatures']
    
    if verbose:
        log(f"\nElbow detection:")
        log(f"  Elbow found at: {elbow_n_sigs} signatures")
        log(f"  Error at elbow: {scree_data[elbow_idx]['frobenius_error']:.2f}")
        log(f"  Relative error: {100*scree_data[elbow_idx]['relative_error']:.2f}%")
        log(f"  Explained variance: {100*scree_data[elbow_idx]['explained_variance']:.2f}%")
    
    # -------------------------------------------------------------------------
    # Step 4: Validate elbow with marginal improvement analysis
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"STEP 4: Validating elbow point (marginal improvement analysis)")
        log(f"{'='*70}")
    
    # Calculate marginal improvement for each step
    marginal_improvements = []
    for i in range(1, len(scree_data)):
        prev_error = scree_data[i-1]['frobenius_error']
        curr_error = scree_data[i]['frobenius_error']
        improvement = 100 * (prev_error - curr_error) / prev_error
        marginal_improvements.append(improvement)
    
    # Find where improvement drops below threshold (e.g., < 1%)
    threshold_pct = 1.0
    threshold_n_sigs = len(scree_data)  # default to all
    for i, imp in enumerate(marginal_improvements):
        if imp < threshold_pct:
            threshold_n_sigs = scree_data[i+1]['n_signatures']
            if verbose:
                log(f"\n  Marginal improvement drops below {threshold_pct}% at {threshold_n_sigs} signatures")
            break
    
    # Use the more conservative of elbow and threshold methods
    final_n_sigs = min(elbow_n_sigs, threshold_n_sigs)
    
    if verbose:
        log(f"\n  Final selection: {final_n_sigs} signatures")
        log(f"    (Elbow method: {elbow_n_sigs}, Threshold method: {threshold_n_sigs})")
    
    selected_signatures = scree_data[final_n_sigs - len(core_available)]['signatures']
    
    # -------------------------------------------------------------------------
    # Step 5: Create scree plot visualization
    # -------------------------------------------------------------------------
    if verbose:
        log(f"\n{'='*70}")
        log(f"Creating scree plot visualization")
        log(f"{'='*70}")
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Frobenius error vs. number of signatures
    ax1 = axes[0, 0]
    ax1.plot(n_sigs_array, errors, 'o-', linewidth=2, markersize=6, color='steelblue')
    ax1.axvline(final_n_sigs, color='red', linestyle='--', linewidth=2, 
                label=f'Selected: {final_n_sigs} sigs')
    ax1.set_xlabel('Number of Signatures', fontsize=11)
    ax1.set_ylabel('Frobenius Error', fontsize=11)
    ax1.set_title('Scree Plot: Reconstruction Error', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Panel 2: Explained variance vs. number of signatures
    ax2 = axes[0, 1]
    explained_vars = np.array([s['explained_variance'] for s in scree_data])
    ax2.plot(n_sigs_array, explained_vars * 100, 'o-', linewidth=2, markersize=6, color='darkgreen')
    ax2.axvline(final_n_sigs, color='red', linestyle='--', linewidth=2,
                label=f'Selected: {final_n_sigs} sigs')
    ax2.axhline(90, color='orange', linestyle=':', alpha=0.5, label='90%')
    ax2.set_xlabel('Number of Signatures', fontsize=11)
    ax2.set_ylabel('Explained Variance (%)', fontsize=11)
    ax2.set_title('Explained Variance', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0, 100])
    
    # Panel 3: Marginal improvement per signature added
    ax3 = axes[1, 0]
    ax3.bar(n_sigs_array[1:], marginal_improvements, color='coral', alpha=0.7, edgecolor='black')
    ax3.axhline(threshold_pct, color='red', linestyle='--', linewidth=2,
                label=f'{threshold_pct}% threshold')
    ax3.axvline(final_n_sigs, color='red', linestyle='--', linewidth=2, alpha=0.5)
    ax3.set_xlabel('Number of Signatures', fontsize=11)
    ax3.set_ylabel('Marginal Improvement (%)', fontsize=11)
    ax3.set_title('Marginal Error Reduction per Signature', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    
    # Panel 4: Second derivative (curvature) for elbow detection
    ax4 = axes[1, 1]
    ax4.plot(n_sigs_array, d2y, 'o-', linewidth=2, markersize=6, color='purple')
    ax4.axhline(0, color='black', linestyle='-', linewidth=0.5)
    ax4.axvline(elbow_n_sigs, color='red', linestyle='--', linewidth=2,
                label=f'Elbow: {elbow_n_sigs} sigs')
    ax4.scatter(elbow_n_sigs, d2y[elbow_idx], color='red', s=150, zorder=5,
                edgecolor='black', linewidth=2)
    ax4.set_xlabel('Number of Signatures', fontsize=11)
    ax4.set_ylabel('Second Derivative (Curvature)', fontsize=11)
    ax4.set_title('Elbow Detection (Maximum Curvature)', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)
    
    plt.suptitle('Signature Selection: Scree Plot Analysis', fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    scree_plot_file = output_path / "signature_selection_scree_plot.png"
    plt.savefig(scree_plot_file, dpi=300, bbox_inches='tight')
    if verbose:
        log(f"Saved scree plot: {scree_plot_file}")
    plt.show()
    
    # -------------------------------------------------------------------------
    # Step 6: Save detailed results
    # -------------------------------------------------------------------------
    results_file = output_path / "signature_selection_scree_analysis.txt"
    
    with open(results_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SIGNATURE SELECTION: SCREE PLOT ANALYSIS\n")
        f.write("="*80 + "\n\n")
        
        f.write("Selection Method: Elbow detection in scree plot\n")
        f.write("Principle: Stop when adding more signatures gives diminishing returns\n\n")
        
        f.write("="*80 + "\n")
        f.write("CANDIDATE RANKING (by individual explanatory power)\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"{'Rank':<6} {'Signature':<12} {'Explained Variance':<20}\n")
        f.write("-"*40 + "\n")
        for i, score in enumerate(candidate_scores, 1):
            f.write(f"{i:<6} {score['signature']:<12} {100*score['explained_variance']:>18.2f}%\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("SCREE PLOT DATA\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"{'N':<4} {'Frobenius Error':<18} {'Rel.Error%':<12} {'Exp.Var%':<12} {'Marginal Imp%':<15}\n")
        f.write("-"*80 + "\n")
        
        for i, data in enumerate(scree_data):
            marg_imp = marginal_improvements[i-1] if i > 0 else 0
            marker = " ← SELECTED" if data['n_signatures'] == final_n_sigs else ""
            f.write(f"{data['n_signatures']:<4} {data['frobenius_error']:<18.2f} "
                   f"{100*data['relative_error']:<11.2f}% {100*data['explained_variance']:<11.2f}% "
                   f"{marg_imp:<14.2f}%{marker}\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("ELBOW DETECTION RESULTS\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Elbow method: {elbow_n_sigs} signatures\n")
        f.write(f"Threshold method (<{threshold_pct}% improvement): {threshold_n_sigs} signatures\n")
        f.write(f"Final selection: {final_n_sigs} signatures (more conservative)\n\n")
        
        f.write("="*80 + "\n")
        f.write("SELECTED SIGNATURES\n")
        f.write("="*80 + "\n\n")
        
        for i, sig in enumerate(selected_signatures, 1):
            marker = "[CORE]" if sig in core_available else "[ADDED]"
            f.write(f"{i:2d}. {sig:<10} {marker}\n")
        
        f.write(f"\nTotal: {len(selected_signatures)} signatures\n")
        f.write(f"Core: {len([s for s in selected_signatures if s in core_available])}\n")
        f.write(f"Added: {len([s for s in selected_signatures if s not in core_available])}\n")
    
    if verbose:
        log(f"Saved analysis: {results_file}")
        
        log(f"\n{'='*70}")
        log(f"SIGNATURE SELECTION COMPLETE")
        log(f"{'='*70}")
        log(f"\nSelected {final_n_sigs} signatures:")
        for i, sig in enumerate(selected_signatures, 1):
            marker = "🔹" if sig in core_available else "✓"
            log(f"  {i:2d}. {marker} {sig}")
        log("")
    
    results = {
        'selected_signatures': selected_signatures,
        'signature_matrix': hnscc_sigs[selected_signatures],
        'n_signatures': final_n_sigs,
        'scree_data': scree_data,
        'candidate_ranking': candidate_scores,
        'elbow_n_sigs': elbow_n_sigs,
        'final_error': scree_data[final_n_sigs - len(core_available)]['frobenius_error'],
        'final_relative_error': scree_data[final_n_sigs - len(core_available)]['relative_error'],
        'final_explained_variance': scree_data[final_n_sigs - len(core_available)]['explained_variance']
    }
    
    return results


# =============================================================================
# STEP 3: SOLVE FOR WEIGHTS USING NNLS
# =============================================================================

def fit_signatures_nnls(mutation_matrix, signature_matrix, verbose=True):
    """
    Fit known signatures to mutation data using Non-Negative Least Squares
    
    Solves: X ≈ H × W
    
    For each cell (column in X), solves:
        min ||X_cell - H × W_cell||²  subject to W_cell ≥ 0
    
    Parameters:
    -----------
    mutation_matrix : pd.DataFrame
        Mutation count matrix (96 contexts × n_cells)
    signature_matrix : pd.DataFrame
        Known signature matrix (96 contexts × k_signatures)
    verbose : bool
        Print progress updates
    
    Returns:
    --------
    dict with:
        - 'weights': pd.DataFrame (k_signatures × n_cells)
        - 'residuals': np.array (residual per cell)
        - 'reconstruction': pd.DataFrame (96 contexts × n_cells)
    """
    if verbose:
        log("="*80)
        log("FITTING SIGNATURES USING NON-NEGATIVE LEAST SQUARES")
        log("="*80)
    
    # Ensure contexts match
    if not all(mutation_matrix.index == signature_matrix.index):
        log("⚠️  Reordering signatures to match mutation matrix contexts...", level="WARNING")
        signature_matrix = signature_matrix.loc[mutation_matrix.index]
    
    n_contexts = mutation_matrix.shape[0]
    n_cells = mutation_matrix.shape[1]
    n_sigs = signature_matrix.shape[1]
    
    if verbose:
        log(f"\nProblem dimensions:")
        log(f"  Mutation contexts: {n_contexts}")
        log(f"  Cells to fit: {n_cells:,}")
        log(f"  Signatures: {n_sigs}")
        log(f"\nSolving: X ({n_contexts}×{n_cells}) ≈ H ({n_contexts}×{n_sigs}) × W ({n_sigs}×{n_cells})")
    
    # Convert to numpy arrays
    X = mutation_matrix.values  # 96 × n_cells
    H = signature_matrix.values  # 96 × n_sigs
    
    # Initialize weight matrix
    W = np.zeros((n_sigs, n_cells))
    residuals = np.zeros(n_cells)
    
    if verbose:
        log(f"\nFitting signatures to {n_cells:,} cells...")
    
    # Solve for each cell
    progress_interval = max(1, n_cells // 20)  # Update every 5%
    
    for i in range(n_cells):
        # NNLS for this cell
        cell_mutations = X[:, i]
        
        try:
            weights, residual = nnls(H, cell_mutations)
            W[:, i] = weights
            residuals[i] = residual
        except Exception as e:
            if verbose and i == 0:
                log(f"⚠️  NNLS failed for cell {i}: {e}", level="WARNING")
            W[:, i] = 0
            residuals[i] = np.inf
        
        # Progress update
        if verbose and (i + 1) % progress_interval == 0:
            pct = 100 * (i + 1) / n_cells
            log(f"  Progress: {i+1:,}/{n_cells:,} cells ({pct:.1f}%)")
    
    if verbose:
        log(f"✓ Fitting complete for all {n_cells:,} cells")
    
    # Create weight DataFrame
    weights_df = pd.DataFrame(
        W,
        index=signature_matrix.columns,
        columns=mutation_matrix.columns
    )
    
    # Calculate reconstruction
    reconstruction = H @ W
    reconstruction_df = pd.DataFrame(
        reconstruction,
        index=mutation_matrix.index,
        columns=mutation_matrix.columns
    )
    
    results = {
        'weights': weights_df,
        'residuals': residuals,
        'reconstruction': reconstruction_df
    }
    
    if verbose:
        log("\n" + "="*80)
        log("SIGNATURE FITTING COMPLETE")
        log("="*80 + "\n")
    
    return results


# =============================================================================
# STEP 4: EVALUATE RECONSTRUCTION QUALITY
# =============================================================================

def evaluate_reconstruction(original_matrix, reconstructed_matrix, verbose=True):
    """
    Evaluate how well H × W approximates X using proper matrix norms
    
    Uses:
    - Frobenius norm for reconstruction error
    - Eckart-Young theorem for optimal low-rank approximation
    - Per-cell correlation metrics
    
    Metrics:
    - Frobenius norm ||X - (H×W)||_F (lower is better)
    - Relative Frobenius error ||X - (H×W)||_F / ||X||_F
    - Eckart-Young optimal rank-k approximation via SVD
    - Per-cell Pearson correlation
    - Per-cell Cosine similarity
    
    Parameters:
    -----------
    original_matrix : pd.DataFrame
        Original mutation matrix X
    reconstructed_matrix : pd.DataFrame
        Reconstructed matrix H × W
    verbose : bool
        Print detailed statistics
    
    Returns:
    --------
    dict with evaluation metrics
    """
    if verbose:
        log("="*80)
        log("EVALUATING RECONSTRUCTION QUALITY")
        log("="*80)
    
    X = original_matrix.values
    X_recon = reconstructed_matrix.values
    
    n_contexts = X.shape[0]
    n_cells = X.shape[1]
    
    if verbose:
        log(f"\nMatrix dimensions: {n_contexts} contexts × {n_cells} cells")
    
    # -------------------------------------------------------------------------
    # Calculate difference matrix (residual)
    # -------------------------------------------------------------------------
    if verbose:
        log("\nCalculating residual matrix (X - X_reconstructed)...")
    
    residual_matrix = X - X_recon
    
    # -------------------------------------------------------------------------
    # Frobenius Norm Analysis
    # -------------------------------------------------------------------------
    if verbose:
        log("\nCalculating Frobenius norms...")
    
    # ||X||_F
    frobenius_norm_original = np.linalg.norm(X, 'fro')
    
    # ||X_recon||_F
    frobenius_norm_reconstructed = np.linalg.norm(X_recon, 'fro')
    
    # ||X - X_recon||_F (reconstruction error)
    frobenius_error = np.linalg.norm(residual_matrix, 'fro')
    
    # Relative error: ||X - X_recon||_F / ||X||_F
    relative_frobenius_error = frobenius_error / frobenius_norm_original if frobenius_norm_original > 0 else 0
    
    if verbose:
        log(f"  ||X||_F (original):         {frobenius_norm_original:.2f}")
        log(f"  ||X_recon||_F:              {frobenius_norm_reconstructed:.2f}")
        log(f"  ||X - X_recon||_F (error):  {frobenius_error:.2f}")
        log(f"  Relative error:             {100*relative_frobenius_error:.2f}%")
    
    # -------------------------------------------------------------------------
    # Eckart-Young Theorem: Optimal Low-Rank Approximation
    # -------------------------------------------------------------------------
    if verbose:
        log("\nApplying Eckart-Young theorem (optimal rank-k approximation via SVD)...")
    
    # Perform SVD on original matrix: X = U Σ V^T
    U, singular_values, Vt = np.linalg.svd(X, full_matrices=False)
    
    # Rank of the matrix
    matrix_rank = np.linalg.matrix_rank(X)
    
    if verbose:
        log(f"  Matrix rank: {matrix_rank}")
        log(f"  Top 10 singular values: {singular_values[:10]}")
    
    # Calculate cumulative explained variance
    total_variance = np.sum(singular_values**2)
    cumulative_variance = np.cumsum(singular_values**2) / total_variance
    
    # Find rank needed for 90%, 95%, 99% reconstruction
    rank_90 = np.searchsorted(cumulative_variance, 0.90) + 1
    rank_95 = np.searchsorted(cumulative_variance, 0.95) + 1
    rank_99 = np.searchsorted(cumulative_variance, 0.99) + 1
    
    if verbose:
        log(f"\n  Eckart-Young optimal ranks:")
        log(f"    90% variance: rank-{rank_90}")
        log(f"    95% variance: rank-{rank_95}")
        log(f"    99% variance: rank-{rank_99}")
    
    # Calculate theoretical minimum error for rank-k approximation
    # More accurate: count non-zero singular values in reconstructed matrix
    _, sv_recon, _ = np.linalg.svd(X_recon, full_matrices=False)
    decomposition_rank_effective = np.sum(sv_recon > 1e-10)
    
    if decomposition_rank_effective < len(singular_values):
        # Theoretical minimum error for this rank
        theoretical_min_error = np.sqrt(np.sum(singular_values[decomposition_rank_effective:]**2))
        
        # How close are we to optimal?
        optimality_ratio = theoretical_min_error / frobenius_error if frobenius_error > 0 else 1.0
        
        if verbose:
            log(f"\n  Effective decomposition rank: {decomposition_rank_effective}")
            log(f"  Theoretical minimum error (Eckart-Young): {theoretical_min_error:.2f}")
            log(f"  Actual error: {frobenius_error:.2f}")
            log(f"  Optimality ratio: {optimality_ratio:.4f}")
            
            if optimality_ratio > 0.95:
                log(f"    ✓ Near-optimal reconstruction (>95% of theoretical best)")
            elif optimality_ratio > 0.80:
                log(f"    ○ Good reconstruction (>80% of theoretical best)")
            else:
                log(f"    ⚠ Suboptimal reconstruction (<80% of theoretical best)")
    else:
        theoretical_min_error = 0
        optimality_ratio = 1.0
        if verbose:
            log(f"\n  Effective decomposition rank: {decomposition_rank_effective} (full rank)")
    
    # -------------------------------------------------------------------------
    # Per-cell metrics
    # -------------------------------------------------------------------------
    if verbose:
        log("\nCalculating per-cell metrics...")
    
    pearson_corrs = []
    spearman_corrs = []
    cosine_sims = []
    cell_frobenius_errors = []
    
    for i in range(n_cells):
        x_cell = X[:, i]
        x_recon_cell = X_recon[:, i]
        
        # Frobenius error for this cell (L2 norm of difference)
        cell_error = np.linalg.norm(x_cell - x_recon_cell)
        cell_frobenius_errors.append(cell_error)
        
        # Pearson correlation
        if x_cell.sum() > 0 and x_recon_cell.sum() > 0:
            r, _ = pearsonr(x_cell, x_recon_cell)
            pearson_corrs.append(r)
            
            # Spearman correlation
            rho, _ = spearmanr(x_cell, x_recon_cell)
            spearman_corrs.append(rho)
            
            # Cosine similarity
            norm_x = np.linalg.norm(x_cell)
            norm_recon = np.linalg.norm(x_recon_cell)
            if norm_x > 0 and norm_recon > 0:
                cos_sim = np.dot(x_cell, x_recon_cell) / (norm_x * norm_recon)
                cosine_sims.append(cos_sim)
            else:
                cosine_sims.append(np.nan)
        else:
            pearson_corrs.append(np.nan)
            spearman_corrs.append(np.nan)
            cosine_sims.append(np.nan)
    
    pearson_corrs = np.array(pearson_corrs)
    spearman_corrs = np.array(spearman_corrs)
    cosine_sims = np.array(cosine_sims)
    cell_frobenius_errors = np.array(cell_frobenius_errors)
    
    # -------------------------------------------------------------------------
    # Additional metrics
    # -------------------------------------------------------------------------
    
    # Reconstruction rate (what % of mutations are explained)
    total_mutations_original = X.sum()
    total_mutations_reconstructed = X_recon.sum()
    reconstruction_rate = total_mutations_reconstructed / total_mutations_original if total_mutations_original > 0 else 0
    
    # Mean Absolute Error per element
    mae = np.mean(np.abs(residual_matrix))
    
    # Root Mean Squared Error per element
    rmse = np.sqrt(np.mean(residual_matrix**2))
    
    # -------------------------------------------------------------------------
    # Summary statistics
    # -------------------------------------------------------------------------
    if verbose:
        log("\n" + "="*60)
        log("RECONSTRUCTION QUALITY METRICS")
        log("="*60)
        
        log("\nMatrix Norm Metrics:")
        log(f"  Frobenius error:          {frobenius_error:.2f}")
        log(f"  Relative Frobenius error: {100*relative_frobenius_error:.2f}%")
        log(f"  MAE (per element):        {mae:.4f}")
        log(f"  RMSE (per element):       {rmse:.4f}")
        
        log("\nEckart-Young Optimal Approximation:")
        log(f"  Theoretical min error:    {theoretical_min_error:.2f}")
        log(f"  Optimality ratio:         {optimality_ratio:.4f}")
        log(f"  Effective rank:           {decomposition_rank_effective}")
        
        log("\nPer-Cell Metrics (averaged across cells):")
        log(f"  Mean cell error (L2):     {np.mean(cell_frobenius_errors):.2f}")
        log(f"  Median cell error:        {np.median(cell_frobenius_errors):.2f}")
        
        log(f"\n  Pearson correlation:")
        log(f"    Mean:   {np.nanmean(pearson_corrs):.4f}")
        log(f"    Median: {np.nanmedian(pearson_corrs):.4f}")
        log(f"    Std:    {np.nanstd(pearson_corrs):.4f}")
        
        log(f"\n  Spearman correlation:")
        log(f"    Mean:   {np.nanmean(spearman_corrs):.4f}")
        log(f"    Median: {np.nanmedian(spearman_corrs):.4f}")
        
        log(f"\n  Cosine similarity:")
        log(f"    Mean:   {np.nanmean(cosine_sims):.4f}")
        log(f"    Median: {np.nanmedian(cosine_sims):.4f}")
        
        log("\nMutation Counts:")
        log(f"  Original total:       {int(total_mutations_original):,}")
        log(f"  Reconstructed total:  {int(total_mutations_reconstructed):,}")
        log(f"  Difference:           {int(total_mutations_original - total_mutations_reconstructed):,}")
        log(f"  Reconstruction rate:  {100*reconstruction_rate:.2f}%")
    
    # Quality assessment based on Frobenius error
    if verbose:
        log("\nQuality Assessment:")
        if relative_frobenius_error < 0.10:
            quality = "EXCELLENT"
        elif relative_frobenius_error < 0.20:
            quality = "GOOD"
        elif relative_frobenius_error < 0.30:
            quality = "MODERATE"
        else:
            quality = "POOR"
        log(f"  Overall quality (Frobenius): {quality}")
        
        # Also assess by correlation
        mean_pearson = np.nanmean(pearson_corrs)
        if mean_pearson > 0.8:
            quality_corr = "EXCELLENT"
        elif mean_pearson > 0.6:
            quality_corr = "GOOD"
        elif mean_pearson > 0.4:
            quality_corr = "MODERATE"
        else:
            quality_corr = "POOR"
        log(f"  Overall quality (Correlation): {quality_corr}")
    
    results = {
        # Matrix difference
        'residual_matrix': residual_matrix,
        
        # Frobenius norms
        'frobenius_norm_original': frobenius_norm_original,
        'frobenius_norm_reconstructed': frobenius_norm_reconstructed,
        'frobenius_error': frobenius_error,
        'relative_frobenius_error': relative_frobenius_error,
        
        # Eckart-Young theorem results
        'singular_values': singular_values,
        'matrix_rank': matrix_rank,
        'decomposition_rank': decomposition_rank_effective,
        'theoretical_min_error': theoretical_min_error,
        'optimality_ratio': optimality_ratio,
        'rank_90pct': rank_90,
        'rank_95pct': rank_95,
        'rank_99pct': rank_99,
        'cumulative_variance': cumulative_variance,
        
        # Per-cell metrics
        'pearson_per_cell': pearson_corrs,
        'spearman_per_cell': spearman_corrs,
        'cosine_per_cell': cosine_sims,
        'cell_frobenius_errors': cell_frobenius_errors,
        
        # Summary statistics
        'mean_pearson': np.nanmean(pearson_corrs),
        'median_pearson': np.nanmedian(pearson_corrs),
        'std_pearson': np.nanstd(pearson_corrs),
        'mean_cosine': np.nanmean(cosine_sims),
        'median_cosine': np.nanmedian(cosine_sims),
        'mean_cell_error': np.mean(cell_frobenius_errors),
        
        # Other metrics
        'mae': mae,
        'rmse': rmse,
        'reconstruction_rate': reconstruction_rate,
        'total_mutations_original': total_mutations_original,
        'total_mutations_reconstructed': total_mutations_reconstructed,
        'quality': quality if verbose else None,
        'quality_correlation': quality_corr if verbose else None
    }
    
    if verbose:
        log("\n" + "="*80)
        log("EVALUATION COMPLETE")
        log("="*80 + "\n")
    
    return results


# =============================================================================
# STEP 5: VISUALIZATION
# =============================================================================

def plot_reconstruction_quality(evaluation_results, output_dir):
    """
    Visualize reconstruction quality metrics including SVD analysis
    """
    log("="*80)
    log("GENERATING RECONSTRUCTION QUALITY PLOTS")
    log("="*80)
    
    fig = plt.figure(figsize=(20, 14))
    
    pearson_corrs = evaluation_results['pearson_per_cell']
    cosine_sims = evaluation_results['cosine_per_cell']
    cell_errors = evaluation_results['cell_frobenius_errors']
    
    # Remove NaN values for plotting
    pearson_valid = pearson_corrs[~np.isnan(pearson_corrs)]
    cosine_valid = cosine_sims[~np.isnan(cosine_sims)]
    
    # -------------------------------------------------------------------------
    # Panel 1: Pearson correlation distribution
    # -------------------------------------------------------------------------
    ax1 = plt.subplot(3, 3, 1)
    ax1.hist(pearson_valid, bins=50, color='steelblue', alpha=0.7, edgecolor='black')
    ax1.axvline(np.nanmean(pearson_corrs), color='red', linestyle='--', 
                linewidth=2, label=f'Mean: {np.nanmean(pearson_corrs):.3f}')
    ax1.axvline(np.nanmedian(pearson_corrs), color='orange', linestyle='--',
                linewidth=2, label=f'Median: {np.nanmedian(pearson_corrs):.3f}')
    ax1.set_xlabel('Pearson Correlation', fontsize=11)
    ax1.set_ylabel('Number of Cells', fontsize=11)
    ax1.set_title('Per-Cell Reconstruction Quality\n(Pearson Correlation)', 
                  fontsize=12, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(axis='y', alpha=0.3)
    
    # -------------------------------------------------------------------------
    # Panel 2: Cosine similarity distribution
    # -------------------------------------------------------------------------
    ax2 = plt.subplot(3, 3, 2)
    ax2.hist(cosine_valid, bins=50, color='darkgreen', alpha=0.7, edgecolor='black')
    ax2.axvline(np.nanmean(cosine_sims), color='red', linestyle='--',
                linewidth=2, label=f'Mean: {np.nanmean(cosine_sims):.3f}')
    ax2.axvline(np.nanmedian(cosine_sims), color='orange', linestyle='--',
                linewidth=2, label=f'Median: {np.nanmedian(cosine_sims):.3f}')
    ax2.set_xlabel('Cosine Similarity', fontsize=11)
    ax2.set_ylabel('Number of Cells', fontsize=11)
    ax2.set_title('Per-Cell Reconstruction Quality\n(Cosine Similarity)',
                  fontsize=12, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(axis='y', alpha=0.3)
    
    # -------------------------------------------------------------------------
    # Panel 3: Per-cell Frobenius error distribution
    # -------------------------------------------------------------------------
    ax3 = plt.subplot(3, 3, 3)
    ax3.hist(cell_errors, bins=50, color='coral', alpha=0.7, edgecolor='black')
    ax3.axvline(np.mean(cell_errors), color='red', linestyle='--',
                linewidth=2, label=f'Mean: {np.mean(cell_errors):.2f}')
    ax3.axvline(np.median(cell_errors), color='orange', linestyle='--',
                linewidth=2, label=f'Median: {np.median(cell_errors):.2f}')
    ax3.set_xlabel('Per-Cell Frobenius Error (L2 norm)', fontsize=11)
    ax3.set_ylabel('Number of Cells', fontsize=11)
    ax3.set_title('Per-Cell Reconstruction Error\n(Frobenius Norm)',
                  fontsize=12, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    
    # -------------------------------------------------------------------------
    # Panel 4: Quality categories
    # -------------------------------------------------------------------------
    ax4 = plt.subplot(3, 3, 4)
    
    # Categorize cells
    excellent = np.sum(pearson_valid > 0.8)
    good = np.sum((pearson_valid > 0.6) & (pearson_valid <= 0.8))
    moderate = np.sum((pearson_valid > 0.4) & (pearson_valid <= 0.6))
    poor = np.sum(pearson_valid <= 0.4)
    
    categories = ['Excellent\n(>0.8)', 'Good\n(0.6-0.8)', 'Moderate\n(0.4-0.6)', 'Poor\n(<0.4)']
    counts = [excellent, good, moderate, poor]
    colors_cat = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c']
    
    bars = ax4.bar(categories, counts, color=colors_cat, alpha=0.7, edgecolor='black')
    ax4.set_ylabel('Number of Cells', fontsize=11)
    ax4.set_title('Reconstruction Quality Categories\n(Pearson Correlation)',
                  fontsize=12, fontweight='bold')
    ax4.grid(axis='y', alpha=0.3)
    
    # Add count labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{count}\n({100*count/len(pearson_valid):.1f}%)',
                ha='center', va='bottom', fontsize=9)
    
    # -------------------------------------------------------------------------
    # Panel 5: Pearson vs Cosine scatter
    # -------------------------------------------------------------------------
    ax5 = plt.subplot(3, 3, 5)
    
    # Match up valid values
    valid_mask = ~np.isnan(pearson_corrs) & ~np.isnan(cosine_sims)
    pearson_for_scatter = pearson_corrs[valid_mask]
    cosine_for_scatter = cosine_sims[valid_mask]
    
    ax5.scatter(pearson_for_scatter, cosine_for_scatter, alpha=0.3, s=10, color='purple')
    ax5.plot([0, 1], [0, 1], 'r--', linewidth=1, alpha=0.5, label='y=x')
    ax5.set_xlabel('Pearson Correlation', fontsize=11)
    ax5.set_ylabel('Cosine Similarity', fontsize=11)
    ax5.set_title('Pearson vs Cosine Similarity',
                  fontsize=12, fontweight='bold')
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.3)
    ax5.set_xlim([0, 1])
    ax5.set_ylim([0, 1])
    
    # -------------------------------------------------------------------------
    # Panel 6: Cumulative distribution
    # -------------------------------------------------------------------------
    ax6 = plt.subplot(3, 3, 6)
    
    sorted_pearson = np.sort(pearson_valid)
    cumulative = np.arange(1, len(sorted_pearson) + 1) / len(sorted_pearson)
    
    ax6.plot(sorted_pearson, cumulative, linewidth=2, color='steelblue')
    ax6.axvline(0.5, color='orange', linestyle='--', alpha=0.5, label='r=0.5')
    ax6.axvline(0.7, color='red', linestyle='--', alpha=0.5, label='r=0.7')
    ax6.set_xlabel('Pearson Correlation', fontsize=11)
    ax6.set_ylabel('Cumulative Fraction of Cells', fontsize=11)
    ax6.set_title('Cumulative Distribution\n(Pearson Correlation)',
                  fontsize=12, fontweight='bold')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)
    ax6.set_xlim([0, 1])
    ax6.set_ylim([0, 1])
    
    # -------------------------------------------------------------------------
    # Panel 7: SVD Scree Plot (Singular Values)
    # -------------------------------------------------------------------------
    ax7 = plt.subplot(3, 3, 7)
    
    singular_values = evaluation_results['singular_values']
    n_show = min(50, len(singular_values))  # Show first 50
    
    ax7.plot(range(1, n_show + 1), singular_values[:n_show], 'o-', 
             linewidth=2, markersize=4, color='navy')
    ax7.set_xlabel('Rank', fontsize=11)
    ax7.set_ylabel('Singular Value', fontsize=11)
    ax7.set_title('SVD Scree Plot\n(Singular Values)',
                  fontsize=12, fontweight='bold')
    ax7.grid(True, alpha=0.3)
    ax7.set_yscale('log')
    
    # Mark decomposition rank
    decomp_rank = evaluation_results['decomposition_rank']
    if decomp_rank <= n_show:
        ax7.axvline(decomp_rank, color='red', linestyle='--', linewidth=2,
                   label=f'Decomp. rank: {decomp_rank}')
        ax7.legend(fontsize=9)
    
    # -------------------------------------------------------------------------
    # Panel 8: Cumulative Variance Explained (Eckart-Young)
    # -------------------------------------------------------------------------
    ax8 = plt.subplot(3, 3, 8)
    
    cumulative_variance = evaluation_results['cumulative_variance']
    
    ax8.plot(range(1, len(cumulative_variance) + 1), cumulative_variance * 100,
             linewidth=2, color='darkgreen')
    ax8.axhline(90, color='orange', linestyle='--', alpha=0.5, label='90%')
    ax8.axhline(95, color='red', linestyle='--', alpha=0.5, label='95%')
    ax8.axhline(99, color='purple', linestyle='--', alpha=0.5, label='99%')
    
    # Mark optimal ranks
    rank_90 = evaluation_results['rank_90pct']
    rank_95 = evaluation_results['rank_95pct']
    rank_99 = evaluation_results['rank_99pct']
    
    ax8.axvline(rank_90, color='orange', linestyle=':', linewidth=1.5, alpha=0.7)
    ax8.axvline(rank_95, color='red', linestyle=':', linewidth=1.5, alpha=0.7)
    ax8.axvline(rank_99, color='purple', linestyle=':', linewidth=1.5, alpha=0.7)
    
    ax8.set_xlabel('Rank', fontsize=11)
    ax8.set_ylabel('Cumulative Variance Explained (%)', fontsize=11)
    ax8.set_title('Eckart-Young Optimal Rank\n(Cumulative Variance)',
                  fontsize=12, fontweight='bold')
    ax8.legend(fontsize=9, loc='lower right')
    ax8.grid(True, alpha=0.3)
    ax8.set_xlim([0, min(100, len(cumulative_variance))])
    ax8.set_ylim([0, 100])
    
    # -------------------------------------------------------------------------
    # Panel 9: Summary statistics
    # -------------------------------------------------------------------------
    ax9 = plt.subplot(3, 3, 9)
    ax9.axis('off')
    
    summary_text = f"""
    RECONSTRUCTION METRICS SUMMARY
    {'='*45}
    
    FROBENIUS NORMS:
      ||X||_F:           {evaluation_results['frobenius_norm_original']:>10.2f}
      ||X_recon||_F:     {evaluation_results['frobenius_norm_reconstructed']:>10.2f}
      ||Error||_F:       {evaluation_results['frobenius_error']:>10.2f}
      Relative error:    {100*evaluation_results['relative_frobenius_error']:>9.2f}%
    
    ECKART-YOUNG THEOREM:
      Matrix rank:       {evaluation_results['matrix_rank']:>10}
      Decomp. rank:      {evaluation_results['decomposition_rank']:>10}
      Theoretical min:   {evaluation_results['theoretical_min_error']:>10.2f}
      Optimality:        {evaluation_results['optimality_ratio']:>10.4f}
    
    CORRELATION METRICS:
      Mean Pearson:      {evaluation_results['mean_pearson']:>10.4f}
      Median Pearson:    {evaluation_results['median_pearson']:>10.4f}
      Mean Cosine:       {evaluation_results['mean_cosine']:>10.4f}
    
    MUTATION COUNTS:
      Original:    {int(evaluation_results['total_mutations_original']):>15,}
      Reconstructed:{int(evaluation_results['total_mutations_reconstructed']):>14,}
      Difference:  {int(evaluation_results['total_mutations_original'] - evaluation_results['total_mutations_reconstructed']):>15,}
    
    QUALITY ASSESSMENT:
      Frobenius:   {evaluation_results['quality']:>15s}
      Correlation: {evaluation_results['quality_correlation']:>15s}
    """
    
    ax9.text(0.05, 0.5, summary_text, fontsize=9, family='monospace',
             verticalalignment='center', transform=ax9.transAxes)
    
    plt.suptitle(f'Reconstruction Quality: H × W ≈ X ({len(pearson_valid):,} cells)',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    save_figure(fig, output_dir, 'reconstruction_quality.png')
    plt.show()
    
    log("\n" + "="*80)
    log("QUALITY PLOTS COMPLETE")
    log("="*80 + "\n")


def plot_signature_weights_summary(weights_df, output_dir, top_n=15):
    """
    Visualize distribution of signature weights across cells
    """
    log("="*80)
    log("GENERATING SIGNATURE WEIGHT PLOTS")
    log("="*80)
    
    n_sigs = weights_df.shape[0]
    n_cells = weights_df.shape[1]
    
    fig = plt.figure(figsize=(18, 10))
    
    # -------------------------------------------------------------------------
    # Panel 1: Mean weights per signature
    # -------------------------------------------------------------------------
    ax1 = plt.subplot(2, 2, 1)
    
    mean_weights = weights_df.mean(axis=1).sort_values(ascending=False)
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(mean_weights)))
    bars = ax1.barh(range(len(mean_weights)), mean_weights.values, color=colors)
    ax1.set_yticks(range(len(mean_weights)))
    ax1.set_yticklabels(mean_weights.index, fontsize=9)
    ax1.set_xlabel('Mean Weight Across Cells', fontsize=11)
    ax1.set_title(f'Average Signature Activity\n({n_cells:,} cells)',
                  fontsize=12, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    ax1.invert_yaxis()
    
    # -------------------------------------------------------------------------
    # Panel 2: Frequency of signature detection
    # -------------------------------------------------------------------------
    ax2 = plt.subplot(2, 2, 2)
    
    # Count cells with weight > 0.1 for each signature
    detection_threshold = 0.1
    detection_freq = (weights_df > detection_threshold).sum(axis=1) / n_cells * 100
    detection_freq = detection_freq.sort_values(ascending=False)
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(detection_freq)))
    bars = ax2.barh(range(len(detection_freq)), detection_freq.values, color=colors)
    ax2.set_yticks(range(len(detection_freq)))
    ax2.set_yticklabels(detection_freq.index, fontsize=9)
    ax2.set_xlabel('% of Cells with Weight > 0.1', fontsize=11)
    ax2.set_title(f'Signature Detection Frequency\n(threshold = {detection_threshold})',
                  fontsize=12, fontweight='bold')
    ax2.grid(axis='x', alpha=0.3)
    ax2.invert_yaxis()
    
    # -------------------------------------------------------------------------
    # Panel 3: Weight distribution heatmap (top signatures)
    # -------------------------------------------------------------------------
    ax3 = plt.subplot(2, 2, 3)
    
    # Select top N signatures by mean weight
    top_sigs = mean_weights.head(min(top_n, n_sigs)).index
    weights_top = weights_df.loc[top_sigs]
    
    # Sample cells if too many
    if n_cells > 1000:
        cell_sample = np.random.choice(weights_df.columns, 1000, replace=False)
        weights_plot = weights_top[cell_sample]
        title_suffix = f" (1000 random cells)"
    else:
        weights_plot = weights_top
        title_suffix = f" ({n_cells} cells)"
    
    im = ax3.imshow(weights_plot.values, aspect='auto', cmap='YlOrRd', 
                    interpolation='nearest')
    ax3.set_yticks(range(len(top_sigs)))
    ax3.set_yticklabels(top_sigs, fontsize=9)
    ax3.set_xlabel('Cells (sampled)', fontsize=11)
    ax3.set_title(f'Signature Weight Heatmap (Top {len(top_sigs)}){title_suffix}',
                  fontsize=12, fontweight='bold')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax3)
    cbar.set_label('Weight', fontsize=10)
    
    # -------------------------------------------------------------------------
    # Panel 4: Total mutations per cell
    # -------------------------------------------------------------------------
    ax4 = plt.subplot(2, 2, 4)
    
    total_mutations_per_cell = weights_df.sum(axis=0)
    
    ax4.hist(total_mutations_per_cell, bins=50, color='teal', 
             alpha=0.7, edgecolor='black')
    ax4.axvline(total_mutations_per_cell.mean(), color='red', 
                linestyle='--', linewidth=2, 
                label=f'Mean: {total_mutations_per_cell.mean():.1f}')
    ax4.axvline(total_mutations_per_cell.median(), color='orange',
                linestyle='--', linewidth=2,
                label=f'Median: {total_mutations_per_cell.median():.1f}')
    ax4.set_xlabel('Total Signature Weight (Sum Across Signatures)', fontsize=11)
    ax4.set_ylabel('Number of Cells', fontsize=11)
    ax4.set_title('Total Signature Activity Per Cell',
                  fontsize=12, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    save_figure(fig, output_dir, 'signature_weights_summary.png')
    plt.show()
    
    log("\n" + "="*80)
    log("WEIGHT PLOTS COMPLETE")
    log("="*80 + "\n")


# =============================================================================
# STEP 6: SAVE RESULTS
# =============================================================================

def save_refitting_results(weights_df, reconstruction_df, evaluation_results,
                           signature_matrix, mutation_matrix, output_dir):
    """
    Save all refitting results to files
    """
    log("="*80)
    log("SAVING REFITTING RESULTS")
    log("="*80)
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # -------------------------------------------------------------------------
    # 1. Signature weights (W matrix)
    # -------------------------------------------------------------------------
    weights_file = output_path / "signature_weights_per_cell.txt"
    weights_df.to_csv(weights_file, sep='\t', float_format='%.6f')
    log(f"Saved weights: {weights_file}")
    
    # -------------------------------------------------------------------------
    # 2. Reconstructed mutation matrix (H × W)
    # -------------------------------------------------------------------------
    recon_file = output_path / "reconstructed_mutation_matrix.txt"
    reconstruction_df.to_csv(recon_file, sep='\t', float_format='%.2f')
    log(f"Saved reconstruction: {recon_file}")
    
    # -------------------------------------------------------------------------
    # 3. Residual matrix (X - X_reconstructed)
    # -------------------------------------------------------------------------
    residual_matrix = evaluation_results['residual_matrix']
    residual_df = pd.DataFrame(
        residual_matrix,
        index=mutation_matrix.index,
        columns=mutation_matrix.columns
    )
    residual_file = output_path / "residual_matrix_X_minus_Xrecon.txt"
    residual_df.to_csv(residual_file, sep='\t', float_format='%.2f')
    log(f"Saved residual matrix: {residual_file}")
    
    # -------------------------------------------------------------------------
    # 4. Evaluation metrics
    # -------------------------------------------------------------------------
    eval_file = output_path / "reconstruction_evaluation.txt"
    
    with open(eval_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SIGNATURE REFITTING EVALUATION\n")
        f.write("="*80 + "\n\n")
        
        f.write("Model:\n")
        f.write(f"  X ≈ H × W\n\n")
        f.write(f"Dimensions:\n")
        f.write(f"  X: {mutation_matrix.shape[0]} contexts × {mutation_matrix.shape[1]} cells\n")
        f.write(f"  H: {signature_matrix.shape[0]} contexts × {signature_matrix.shape[1]} signatures (FIXED)\n")
        f.write(f"  W: {weights_df.shape[0]} signatures × {weights_df.shape[1]} cells (FITTED)\n\n")
        
        f.write("="*80 + "\n")
        f.write("MATRIX NORM METRICS\n")
        f.write("="*80 + "\n\n")
        
        f.write("Frobenius Norms:\n")
        f.write(f"  ||X||_F (original):           {evaluation_results['frobenius_norm_original']:.2f}\n")
        f.write(f"  ||X_recon||_F:                {evaluation_results['frobenius_norm_reconstructed']:.2f}\n")
        f.write(f"  ||X - X_recon||_F (error):    {evaluation_results['frobenius_error']:.2f}\n")
        f.write(f"  Relative error:               {100*evaluation_results['relative_frobenius_error']:.2f}%\n\n")
        
        f.write("Element-wise Metrics:\n")
        f.write(f"  MAE (per element):            {evaluation_results['mae']:.4f}\n")
        f.write(f"  RMSE (per element):           {evaluation_results['rmse']:.4f}\n\n")
        
        f.write("="*80 + "\n")
        f.write("ECKART-YOUNG THEOREM (OPTIMAL LOW-RANK APPROXIMATION)\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Matrix Properties:\n")
        f.write(f"  Full matrix rank:             {evaluation_results['matrix_rank']}\n")
        f.write(f"  Effective decomposition rank: {evaluation_results['decomposition_rank']}\n\n")
        
        f.write(f"Optimal Ranks for Variance Thresholds:\n")
        f.write(f"  90% variance: rank-{evaluation_results['rank_90pct']}\n")
        f.write(f"  95% variance: rank-{evaluation_results['rank_95pct']}\n")
        f.write(f"  99% variance: rank-{evaluation_results['rank_99pct']}\n\n")
        
        f.write(f"Reconstruction Optimality:\n")
        f.write(f"  Theoretical min error:        {evaluation_results['theoretical_min_error']:.2f}\n")
        f.write(f"  Actual error:                 {evaluation_results['frobenius_error']:.2f}\n")
        f.write(f"  Optimality ratio:             {evaluation_results['optimality_ratio']:.4f}\n")
        
        if evaluation_results['optimality_ratio'] > 0.95:
            f.write(f"  Assessment: Near-optimal (>95% of theoretical best)\n\n")
        elif evaluation_results['optimality_ratio'] > 0.80:
            f.write(f"  Assessment: Good (>80% of theoretical best)\n\n")
        else:
            f.write(f"  Assessment: Suboptimal (<80% of theoretical best)\n\n")
        
        f.write("="*80 + "\n")
        f.write("PER-CELL QUALITY METRICS\n")
        f.write("="*80 + "\n\n")
        
        f.write("Per-Cell Metrics (averaged):\n")
        f.write(f"  Mean cell error (L2 norm):    {evaluation_results['mean_cell_error']:.2f}\n\n")
        
        f.write(f"  Mean Pearson correlation:     {evaluation_results['mean_pearson']:.4f}\n")
        f.write(f"  Median Pearson correlation:   {evaluation_results['median_pearson']:.4f}\n")
        f.write(f"  Std Pearson correlation:      {evaluation_results['std_pearson']:.4f}\n\n")
        
        f.write(f"  Mean Cosine similarity:       {evaluation_results['mean_cosine']:.4f}\n")
        f.write(f"  Median Cosine similarity:     {evaluation_results['median_cosine']:.4f}\n\n")
        
        f.write("="*80 + "\n")
        f.write("MUTATION COUNTS\n")
        f.write("="*80 + "\n\n")
        
        f.write("Mutation Counts:\n")
        f.write(f"  Original total:       {int(evaluation_results['total_mutations_original']):>12,}\n")
        f.write(f"  Reconstructed total:  {int(evaluation_results['total_mutations_reconstructed']):>12,}\n")
        f.write(f"  Difference:           {int(evaluation_results['total_mutations_original'] - evaluation_results['total_mutations_reconstructed']):>12,}\n")
        f.write(f"  Reconstruction rate:  {100*evaluation_results['reconstruction_rate']:>11.2f}%\n\n")
        
        f.write("="*80 + "\n")
        f.write("OVERALL QUALITY ASSESSMENT\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"Based on Frobenius norm:   {evaluation_results['quality']}\n")
        f.write(f"Based on correlation:      {evaluation_results['quality_correlation']}\n")
    
    log(f"Saved evaluation: {eval_file}")
    
    # -------------------------------------------------------------------------
    # 5. SVD analysis (singular values and variance explained)
    # -------------------------------------------------------------------------
    svd_file = output_path / "svd_analysis_singular_values.txt"
    
    singular_values = evaluation_results['singular_values']
    cumulative_variance = evaluation_results['cumulative_variance']
    
    svd_df = pd.DataFrame({
        'rank': np.arange(1, len(singular_values) + 1),
        'singular_value': singular_values,
        'variance_explained': (singular_values**2) / np.sum(singular_values**2),
        'cumulative_variance': cumulative_variance
    })
    
    svd_df.to_csv(svd_file, sep='\t', index=False, float_format='%.6f')
    log(f"Saved SVD analysis: {svd_file}")
    
    # -------------------------------------------------------------------------
    # 6. Per-cell quality metrics
    # -------------------------------------------------------------------------
    quality_df = pd.DataFrame({
        'cell_id': mutation_matrix.columns,
        'frobenius_error': evaluation_results['cell_frobenius_errors'],
        'pearson_correlation': evaluation_results['pearson_per_cell'],
        'cosine_similarity': evaluation_results['cosine_per_cell']
    })
    
    quality_file = output_path / "per_cell_quality_metrics.txt"
    quality_df.to_csv(quality_file, sep='\t', index=False, float_format='%.4f')
    log(f"Saved per-cell metrics: {quality_file}")
    
    # -------------------------------------------------------------------------
    # 7. Signature weight summary
    # -------------------------------------------------------------------------
    mean_weights = weights_df.mean(axis=1).sort_values(ascending=False)
    median_weights = weights_df.median(axis=1).sort_values(ascending=False)
    max_weights = weights_df.max(axis=1).sort_values(ascending=False)
    detection_freq = (weights_df > 0.1).sum(axis=1) / weights_df.shape[1] * 100
    
    summary_df = pd.DataFrame({
        'signature': mean_weights.index,
        'mean_weight': mean_weights.values,
        'median_weight': median_weights.loc[mean_weights.index].values,
        'max_weight': max_weights.loc[mean_weights.index].values,
        'detection_frequency_pct': detection_freq.loc[mean_weights.index].values
    })
    
    summary_file = output_path / "signature_weight_summary.txt"
    summary_df.to_csv(summary_file, sep='\t', index=False, float_format='%.4f')
    log(f"Saved weight summary: {summary_file}")
    
    log("\n" + "="*80)
    log("ALL RESULTS SAVED")
    log("="*80 + "\n")


# =============================================================================
# MAIN PIPELINE
# =============================================================================

def run_signature_refitting_pipeline(mutation_matrix_file, cosmic_file, output_dir,
                                     mutation_threshold=0, use_scree_plot=True,
                                     core_signatures=None, candidate_order=None):
    """
    Complete semi-supervised signature refitting pipeline
    
    Steps:
    1. Extract HNSCC-relevant COSMIC signatures (H matrix)
    2. Load and optionally filter mutation matrix by mutation count
    3. Select optimal signatures using scree plot elbow detection
    4. Fit signatures to mutation data using NNLS (solve for W)
    5. Evaluate reconstruction quality (Frobenius norm, Eckart-Young)
    6. Visualize results
    7. Save everything
    
    Parameters:
    -----------
    mutation_matrix_file : str
        Path to mutation count matrix (96 contexts × cells)
    cosmic_file : str
        Path to COSMIC signature database
    output_dir : str
        Output directory for all results
    mutation_threshold : int
        Minimum mutations per cell to include (default: 0, no filtering)
        Set to >0 to filter out sparse cells
    use_scree_plot : bool
        If True, use scree plot for signature selection (default: True)
        If False, use all HNSCC signatures
    core_signatures : list, optional
        Core signatures to always include (default: ['SBS2', 'SBS3', 'SBS5'])
    candidate_order : list, optional
        Order to try adding candidate signatures (by mean weight)
    
    Returns:
    --------
    dict with all results
    """
    
    print("\n" + "="*80)
    print("SEMI-SUPERVISED SIGNATURE REFITTING PIPELINE")
    print("Using Known HNSCC COSMIC Signatures")
    print("="*80)
    print(f"\nMutation matrix: {mutation_matrix_file}")
    print(f"COSMIC database: {cosmic_file}")
    print(f"Output: {output_dir}")
    print(f"Mutation threshold: {'No filtering' if mutation_threshold == 0 else f'≥{mutation_threshold} mutations/cell'}")
    print(f"Scree plot selection: {'Enabled' if use_scree_plot else 'Disabled (using all signatures)'}")
    print("="*80 + "\n")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # -------------------------------------------------------------------------
    # STEP 1: Extract HNSCC signatures
    # -------------------------------------------------------------------------
    hnscc_sigs = extract_hnscc_signatures(cosmic_file, output_dir)
    
    # -------------------------------------------------------------------------
    # STEP 2: Load mutation matrix
    # -------------------------------------------------------------------------
    log("="*80)
    log("LOADING MUTATION MATRIX")
    log("="*80)
    
    log(f"\nLoading: {mutation_matrix_file}")
    mutation_matrix_full = pd.read_csv(mutation_matrix_file, sep='\t', index_col=0)
    
    log(f"Original matrix: {mutation_matrix_full.shape[0]} contexts × {mutation_matrix_full.shape[1]:,} cells")
    
    # -------------------------------------------------------------------------
    # STEP 2B: Filter by mutation count (if threshold > 0)
    # -------------------------------------------------------------------------
    if mutation_threshold > 0:
        log(f"\n{'='*70}")
        log(f"FILTERING CELLS BY MUTATION COUNT (threshold ≥{mutation_threshold})")
        log(f"{'='*70}")
        
        # Calculate mutations per cell
        mutations_per_cell = mutation_matrix_full.sum(axis=0)
        
        log(f"\nOriginal statistics:")
        log(f"  Total cells: {len(mutations_per_cell):,}")
        log(f"  Mean mutations/cell: {mutations_per_cell.mean():.2f}")
        log(f"  Median mutations/cell: {mutations_per_cell.median():.0f}")
        log(f"  Min: {mutations_per_cell.min():.0f}, Max: {mutations_per_cell.max():.0f}")
        
        # Filter
        cells_keep = mutations_per_cell >= mutation_threshold
        mutation_matrix = mutation_matrix_full.loc[:, cells_keep]
        mutations_filtered = mutations_per_cell[cells_keep]
        
        n_kept = mutation_matrix.shape[1]
        n_removed = mutation_matrix_full.shape[1] - n_kept
        pct_kept = 100 * n_kept / mutation_matrix_full.shape[1]
        
        log(f"\nFiltered statistics:")
        log(f"  Cells kept: {n_kept:,} ({pct_kept:.2f}%)")
        log(f"  Cells removed: {n_removed:,}")
        log(f"  Mean mutations/cell: {mutations_filtered.mean():.2f}")
        log(f"  Median mutations/cell: {mutations_filtered.median():.0f}")
        
        # Calculate sparsity
        sparsity = (mutation_matrix == 0).sum().sum() / (mutation_matrix.shape[0] * mutation_matrix.shape[1])
        log(f"  Sparsity: {100*sparsity:.2f}%")
        
        # Save filtered matrix
        filtered_matrix_file = output_path / f"filtered_mutation_matrix_min{mutation_threshold}muts.txt"
        mutation_matrix.to_csv(filtered_matrix_file, sep='\t', float_format='%.0f')
        log(f"\nSaved filtered matrix: {filtered_matrix_file}")
        
        log("\n" + "="*80)
        log("FILTERING COMPLETE")
        log("="*80 + "\n")
    else:
        log(f"\nNo filtering applied (using all {mutation_matrix_full.shape[1]:,} cells)")
        mutation_matrix = mutation_matrix_full
        mutations_per_cell = mutation_matrix.sum(axis=0)
        log(f"Mean mutations/cell: {mutations_per_cell.mean():.2f}")
        log(f"Median mutations/cell: {mutations_per_cell.median():.0f}")
    
    log(f"\nTotal mutations: {int(mutation_matrix.sum().sum()):,}")
    
    log("\n" + "="*80)
    log("MUTATION MATRIX LOADED")
    log("="*80 + "\n")
    
    # -------------------------------------------------------------------------
    # STEP 3: Signature selection via scree plot
    # -------------------------------------------------------------------------
    if use_scree_plot and (core_signatures is not None or candidate_order is not None):
        # Set defaults if only partially provided
        if core_signatures is None:
            core_signatures = ['SBS5']  # Minimal default - just SBS5
        
        if candidate_order is None:
            # Use all remaining HNSCC signatures as candidates
            candidate_order = [sig for sig in hnscc_sigs.columns if sig not in core_signatures]
        
        # Handle SBS40 variants - add SBS40a if it exists and SBS40 was requested
        if 'SBS40a' in hnscc_sigs.columns:
            if 'SBS40' in core_signatures:
                core_signatures = [s if s != 'SBS40' else 'SBS40a' for s in core_signatures]
            if 'SBS40' in candidate_order:
                candidate_order = [s if s != 'SBS40' else 'SBS40a' for s in candidate_order]
            elif 'SBS40a' not in candidate_order and 'SBS40a' not in core_signatures:
                candidate_order = ['SBS40a'] + candidate_order
        
        # Run scree plot selection
        selection_results = select_signatures_via_scree_plot(
            mutation_matrix,
            hnscc_sigs,
            core_signatures,
            candidate_order,
            output_dir,
            max_signatures=15,
            verbose=True
        )
        
        # Use selected signatures
        hnscc_sigs_final = selection_results['signature_matrix']
        
        log(f"Using {selection_results['n_signatures']} selected signatures via scree plot")
    else:
        # Use all HNSCC signatures (default behavior)
        hnscc_sigs_final = hnscc_sigs
        selection_results = None
        log(f"Using all {hnscc_sigs.shape[1]} HNSCC signatures (no scree plot selection)")
    
    # -------------------------------------------------------------------------
    # STEP 4: Fit signatures using NNLS
    # -------------------------------------------------------------------------
    fitting_results = fit_signatures_nnls(mutation_matrix, hnscc_sigs_final, verbose=True)
    
    weights_df = fitting_results['weights']
    reconstruction_df = fitting_results['reconstruction']
    
    # -------------------------------------------------------------------------
    # STEP 5: Evaluate reconstruction quality
    # -------------------------------------------------------------------------
    evaluation_results = evaluate_reconstruction(
        mutation_matrix, 
        reconstruction_df,
        verbose=True
    )
    
    # -------------------------------------------------------------------------
    # STEP 6: Visualize results
    # -------------------------------------------------------------------------
    plot_reconstruction_quality(evaluation_results, output_dir)
    plot_signature_weights_summary(weights_df, output_dir)
    
    # -------------------------------------------------------------------------
    # STEP 7: Save results
    # -------------------------------------------------------------------------
    save_refitting_results(
        weights_df,
        reconstruction_df,
        evaluation_results,
        hnscc_sigs_final,
        mutation_matrix,
        output_dir
    )
    
    # -------------------------------------------------------------------------
    # Final summary
    # -------------------------------------------------------------------------
    print("\n" + "="*80)
    print("PIPELINE COMPLETE!")
    print("="*80)
    
    print(f"\nSignatures fitted: {hnscc_sigs_final.shape[1]}")
    if use_scree_plot and selection_results:
        print(f"  Core signatures: {len(core_signatures)}")
        print(f"  Selected via scree plot: {', '.join(selection_results['selected_signatures'])}")
    
    print(f"\nCells analyzed: {mutation_matrix.shape[1]:,}")
    if mutation_threshold > 0:
        print(f"  Original cells: {mutation_matrix_full.shape[1]:,}")
        print(f"  Filtered (≥{mutation_threshold} muts): {mutation_matrix.shape[1]:,}")
    
    print(f"\nReconstruction Quality:")
    print(f"  Mean Pearson correlation: {evaluation_results['mean_pearson']:.4f}")
    print(f"  Frobenius error: {evaluation_results['frobenius_error']:.2f}")
    print(f"  Relative error: {100*evaluation_results['relative_frobenius_error']:.2f}%")
    print(f"  Optimality ratio: {evaluation_results['optimality_ratio']:.4f}")
    print(f"  Overall quality (Frobenius): {evaluation_results['quality']}")
    print(f"  Overall quality (Correlation): {evaluation_results['quality_correlation']}")
    
    print(f"\nTop 5 signatures by mean weight:")
    mean_weights = weights_df.mean(axis=1).sort_values(ascending=False)
    for i, (sig, weight) in enumerate(mean_weights.head(5).items(), 1):
        print(f"  {i}. {sig}: {weight:.4f}")
    
    print(f"\nAll results saved to: {output_dir}")
    print("="*80 + "\n")
    
    return {
        'hnscc_signatures': hnscc_sigs_final,
        'mutation_matrix': mutation_matrix,
        'weights': weights_df,
        'reconstruction': reconstruction_df,
        'evaluation': evaluation_results,
        'selection': selection_results,
        'mutation_threshold': mutation_threshold
    }

#%%
# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    
    # =========================================================================
    # CONFIGURATION
    # =========================================================================
    # These should all be defined by the user as inputs for the ClusterCatcher pipeline and read in from the config file
    MUTATION_MATRIX_FILE = os.path.join(working_dir, "SNP_matrix_for_SigProfiler.txt")
    COSMIC_FILE = os.path.join(working_dir, "reference", "COSMIC_v3.4_SBS_GRCh38.txt")
    OUTPUT_DIR = os.path.join(working_dir, "signature_refitting")
    
    # =========================================================================
    # PARAMETERS
    # =========================================================================
    # For the ClusterCatcher pipeline all of these should be set by the user and add in the create_config command with approriate paths.
    # MUTATION FILTERING
    # Set to 0 to use all cells (no filtering)
    # Set to >0 to filter cells with fewer mutations
    MUTATION_THRESHOLD = 0  # Start with no filtering
    
    # SCREE PLOT SELECTION
    # If True: Use scree plot to find optimal number of signatures
    # If False: Use all signatures
    USE_SCREE_PLOT = False  # Recommended: True

    # Core signatures. This should be something that can be set by the user otherwise start with SBS5
    CORE_SIGNATURES = ['SBS2', 'SBS13', 'SBS5']
    
    # Candidate order (by mean weight in dataset)
    CANDIDATE_ORDER = ['SBS1', 'SBS2', 'SBS3', 'SBS5', 'SBS8', 'SBS9',
                       'SBS18', 'SBS40', 'SBS17a', 'SBS17b', 'SBS37', 'SBS41']
    
    # =========================================================================
    # RUN PIPELINE
    # =========================================================================
    
    results = run_signature_refitting_pipeline(
        mutation_matrix_file=MUTATION_MATRIX_FILE,
        cosmic_file=COSMIC_FILE,
        output_dir=OUTPUT_DIR,
        mutation_threshold=MUTATION_THRESHOLD,
        use_scree_plot=USE_SCREE_PLOT,
        core_signatures=CORE_SIGNATURES,
        candidate_order=CANDIDATE_ORDER
    )
    
    print("\n✓ Signature refitting complete!")
    print(f"Check results in: {OUTPUT_DIR}")

#%% Add SBS COSMIC signature weights
import pandas as pd

# Load the signature weights file
sbs_weights = pd.read_csv(os.path.join(working_dir, 'signature_refitting/signature_weights_per_cell.txt'), 
                          sep='\t', index_col=0)

# Each row is a signature, each column is a cell barcode
# Transpose so columns become rows (easier to iterate)
sbs_weights_t = sbs_weights.T

# Add each signature as a new column in adata_pp.obs
for signature in sbs_weights.index:
    # Map signature weights to cells, filling missing with 0
    adata_pp.obs[signature] = adata_pp.obs.index.map(sbs_weights.loc[signature]).fillna(0)

#%% Save your work
for target in list(adata_pp.obs.columns):
    adata_pp.obs[target] = [str(element) for element in adata_pp.obs[target]]

adata_pp.write("adata_final.h5ad")

print(f"Added {len(sbs_weights.index)} SBS signatures to adata_pp.obs")
print(f"Signature names: {list(sbs_weights.index)}")
#%% Generate UMAPs for all SBS signatures and save to figures directory
import os
import matplotlib.pyplot as plt
import scanpy as sc

# Get all SBS signature columns that were added to adata_pp
sbs_signature_columns = [col for col in adata_pp.obs.columns if col.startswith('SBS')]

print(f"Generating UMAPs for {len(sbs_signature_columns)} SBS signatures...")

# Set consistent figure parameters
sc.set_figure_params(scanpy=True, fontsize=25)

# Create a subdirectory for signature UMAPs
signature_umap_dir = os.path.join(figures_dir, 'SBS_signature_UMAPs')
os.makedirs(signature_umap_dir, exist_ok=True)

# Iterate through all SBS signatures and save UMAPs
for signature in sbs_signature_columns:
    # Convert to numeric if stored as string
    if adata_pp.obs[signature].dtype == 'object':
        adata_pp.obs[signature] = pd.to_numeric(adata_pp.obs[signature], errors='coerce').fillna(0)
    
    # Calculate vmax as 95th percentile to handle outliers, minimum of 1 to avoid empty plots
    signature_values = adata_pp.obs[signature].values
    vmax = max(np.percentile(signature_values[signature_values > 0], 95) if (signature_values > 0).any() else 1, 1)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    sc.pl.umap(
        adata_pp, 
        color=signature, 
        size=5, 
        frameon=False, 
        cmap='plasma', 
        ax=ax, 
        show=False,
        vmax=vmax,
        title=f'{signature} Weight'
    )
    
    plt.tight_layout()
    
    # Save figure
    umap_path = os.path.join(signature_umap_dir, f'UMAP_{signature}_weight.png')
    fig.savefig(umap_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"  Saved: {umap_path}")

print(f"\nAll {len(sbs_signature_columns)} SBS signature UMAPs saved to: {signature_umap_dir}")

# Also create a combined summary figure with top signatures by mean weight
top_n = min(9, len(sbs_signature_columns))  # 3x3 grid maximum
if top_n > 0:
    # Calculate mean weights and get top signatures
    mean_weights = {sig: adata_pp.obs[sig].astype(float).mean() for sig in sbs_signature_columns}
    top_signatures = sorted(mean_weights.keys(), key=lambda x: mean_weights[x], reverse=True)[:top_n]
    
    # Determine grid size
    n_cols = 3
    n_rows = (top_n + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))
    axes = axes.flatten() if top_n > 1 else [axes]
    
    for idx, signature in enumerate(top_signatures):
        ax = axes[idx]
        
        # Ensure numeric
        if adata_pp.obs[signature].dtype == 'object':
            adata_pp.obs[signature] = pd.to_numeric(adata_pp.obs[signature], errors='coerce').fillna(0)
        
        signature_values = adata_pp.obs[signature].values
        vmax = max(np.percentile(signature_values[signature_values > 0], 95) if (signature_values > 0).any() else 1, 1)
        
        sc.pl.umap(
            adata_pp, 
            color=signature, 
            size=3, 
            frameon=False, 
            cmap='plasma', 
            ax=ax, 
            show=False,
            vmax=vmax,
            title=f'{signature}\n(mean={mean_weights[signature]:.3f})'
        )
    
    # Hide empty subplots
    for idx in range(top_n, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('Top SBS Signatures by Mean Weight', fontsize=28, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    summary_path = os.path.join(figures_dir, 'SBS_signatures_summary_UMAPs.png')
    fig.savefig(summary_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved summary figure: {summary_path}")

#%% End of the pipeline
import sys
sys.exit(0)
