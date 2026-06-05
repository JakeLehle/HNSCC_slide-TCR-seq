#!/usr/bin/env python

#### Comprehensive SComatic Pipeline for Single-Cell Mutation Analysis ####
#### Optimized for Robustness, Diagnostics, and Performance ####

#%%
import os
import yaml
import pandas as pd
import scanpy as sc
from tqdm import tqdm
import pysam
import subprocess
import shutil
import glob
import pickle
import sys
import threading
import time
import psutil
import numpy as np
from multiprocessing.pool import ThreadPool
from math import ceil
from time import perf_counter
from os import cpu_count

# # Set higher file handle limit for parallel processing
# try:
#     import resource
#     resource.setrlimit(resource.RLIMIT_NOFILE, (65536, 65536))
# except ImportError:
#     print("Resource module not available on this platform")


def create_filtered_bam(input_bam, output_bam, valid_barcodes):
    """
    Create filtered BAM file containing only specified barcodes
    Handles both raw barcodes (from adata) and CB tags with suffix (e.g., '-1')
    """
    print(f"Creating filtered BAM: {output_bam}")
    print(f"Valid barcodes count: {len(valid_barcodes)}")
    
    kept_reads = 0
    total_reads = 0
    
    try:
        with pysam.AlignmentFile(input_bam, "rb") as inbam:
            with pysam.AlignmentFile(output_bam, "wb", template=inbam) as outbam:
                for read in tqdm(inbam, desc="Filtering BAM", miniters=1000):
                    total_reads += 1
                    
                    # Get CB tag and handle different formats
                    if not read.has_tag('CB'):
                        continue
                    
                    cb_tag = read.get_tag('CB')
                    
                    # Extract the raw barcode (before any suffix)
                    raw_barcode = cb_tag.split('-')[0]
                    
                    if raw_barcode in valid_barcodes:
                        outbam.write(read)
                        kept_reads += 1
        
        print(f"Total reads processed: {total_reads:,}")
        print(f"Reads kept: {kept_reads:,}")
        print(f"Reads filtered out: {total_reads - kept_reads:,}")
        
        if kept_reads == 0:
            print("\nWARNING: No reads were kept in the filtered BAM!")
            return False
        else:
            print("Indexing filtered BAM...")
            pysam.index(output_bam)
            return True
            
    except Exception as e:
        print(f"CRITICAL ERROR in create_filtered_bam: {str(e)}")
        return False

def run_scomatic(args):
    """Run SComatic SplitBam for a single sample with enhanced diagnostics"""
    filtered_bam, annotation_file, run_acc, scomatic_dir, scomatic_scripts_dir = args
    cmd = [
        "python",
        os.path.join(scomatic_scripts_dir, "SplitBam", "SplitBamCellTypes.py"),
        "--bam", filtered_bam,
        "--meta", annotation_file,
        "--id", run_acc,
        "--max_nM", "5",
        "--max_NH", "1",
        "--min_MQ", "255",
        "--outdir", scomatic_dir
    ]
    try:
        print(f"Running SComatic for {run_acc}...")
        result = subprocess.run(
            cmd, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True
        )
        print(f"Completed SComatic for {run_acc}")
        # Log first 10 lines of stdout/stderr
        print("STDOUT (head):", '\n'.join(result.stdout.split('\n')[:10]))
        print("STDERR (head):", '\n'.join(result.stderr.split('\n')[:10]))
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: SComatic failed for {run_acc}")
        print(f"Exit code: {e.returncode}")
        print(f"Error message: {e.stderr[:500]}...")  # First 500 chars
        return False
    except Exception as e:
        print(f"UNEXPECTED ERROR in run_scomatic: {str(e)}")
        return False

def run_basecell_counter(args):
    """Run BaseCellCounter with comprehensive validation and diagnostics"""
    bam_file, ref_genome, output_dir, scomatic_scripts_dir = args
    cell_type = os.path.basename(bam_file).split('.')[-2]
    temp_dir = os.path.join(output_dir, f'temp_{cell_type}')
    
    try:
        # Validate inputs before processing
        if not os.path.exists(bam_file) or os.path.getsize(bam_file) == 0:
            print(f"ERROR: BAM file missing or empty: {bam_file}")
            return False
            
        if not os.path.exists(ref_genome) or os.path.getsize(ref_genome) == 0:
            print(f"ERROR: Reference genome missing or empty: {ref_genome}")
            return False
        
        os.makedirs(temp_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "BaseCellCounter", "BaseCellCounter.py"),
            "--bam", bam_file,
            "--ref", ref_genome,
            "--chrom", "all",
            "--out_folder", output_dir,
            "--min_bq", "30",
            "--tmp_dir", temp_dir,
            "--nprocs", "1"
        ]
        
        print(f"\nRunning BaseCellCounter for {cell_type}...")
        print(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Check for output files
        output_files = glob.glob(os.path.join(output_dir, "*.tsv"))
        if not output_files:
            print(f"ERROR: No output files created for {cell_type}")
            return False
            
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: BaseCellCounter failed for {cell_type}")
        print(f"Exit code: {e.returncode}")
        print(f"Error output:\n{e.stderr[:500]}...")
        return False
    except subprocess.TimeoutExpired:
        print(f"\nERROR: BaseCellCounter timed out for {cell_type}")
        return False
    except Exception as e:
        print(f"\nERROR: Unexpected error processing {cell_type}: {str(e)}")
        return False
    finally:
        # Keep temp dirs for debugging failures
        if 'result' in locals() and result.returncode != 0 and os.path.exists(temp_dir):
            print(f"Keeping temp dir for debugging: {temp_dir}")

def run_basecell_counter_parallel(basecell_counter_args, n_workers):
    """Run BaseCellCounter in parallel with comprehensive reporting"""
    if not basecell_counter_args:
        print("\nWARNING: No BAM files to process in BaseCellCounter")
        return []
    
    print("\n" + "="*80)
    print("STEP 2: Running BaseCellCounter")
    print(f"Processing {len(basecell_counter_args)} BAM files with {n_workers} workers")
    print("="*80)
    
    start = perf_counter()
    
    # Validate all BAM files exist before starting
    missing_bams = [bam for bam, *_ in basecell_counter_args if not os.path.exists(bam)]
    if missing_bams:
        print("\nERROR: Missing BAM files:")
        for bam in missing_bams[:5]:
            print(f"- {bam}")
        if len(missing_bams) > 5:
            print(f"... and {len(missing_bams)-5} more")
        raise FileNotFoundError(f"{len(missing_bams)} BAM files not found")
    
    # Run processing
    with ThreadPool(n_workers) as pool:
        chunksize = ceil(len(basecell_counter_args) / (n_workers))
        results = list(pool.map(run_basecell_counter, basecell_counter_args, chunksize=chunksize))
    
    end = perf_counter()
    
    # Generate summary report
    success_count = sum(results)
    failure_count = len(results) - success_count
    duration = end - start
    
    print("\n" + "="*80)
    print("BASE CELL COUNTER SUMMARY:")
    print(f"Total BAMs processed: {len(results)}")
    print(f"Successfully processed: {success_count}")
    print(f"Failures: {failure_count}")
    print(f"Time elapsed: {duration:.2f} seconds")
    print(f"Average time per BAM: {duration/len(results):.2f} seconds")
    
    # Save detailed results
    results_file = os.path.join(os.getcwd(), "basecell_counter_results.csv")
    with open(results_file, 'w') as f:
        f.write("BAM File,Success\n")
        for (bam, *_), success in zip(basecell_counter_args, results):
            f.write(f"{bam},{success}\n")
    print(f"Detailed results saved to: {results_file}")
    print("="*80)
    
    return results

def run_merge_counts(tsv_folder, output_file, scomatic_scripts_dir, sample_id):
    """Run MergeBaseCellCounts.py with validation"""
    try:
        # Check if there are any TSV files to merge
        tsv_files = glob.glob(os.path.join(tsv_folder, "*.tsv"))
        if not tsv_files:
            print(f"WARNING: No TSV files found in {tsv_folder}")
            return False
            
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "MergeCounts", "MergeBaseCellCounts.py"),
            "--tsv_folder", tsv_folder,
            "--outfile", output_file
        ]
        print(f"Merging {len(tsv_files)} count files from {tsv_folder}...")
        
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Verify output was created
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            print(f"ERROR: Merged file not created or empty: {output_file}")
            return False
            
        print(f"Merged counts saved to {output_file}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Failed to merge counts: {e}")
        print(f"Error output: {e.stderr[:500]}...")
        return False
    except Exception as e:
        print(f"ERROR: Unexpected error in merge_counts: {str(e)}")
        return False

def run_variant_calling_step1(args):
    """Run BaseCellCalling step 1 with enhanced validation"""
    merged_counts, output_prefix, ref_genome, scomatic_scripts_dir = args
    step1_output = f"{output_prefix}.calling.step1.tsv"
    
    try:
        # Verify input file
        if not os.path.exists(merged_counts) or os.path.getsize(merged_counts) == 0:
            raise ValueError(f"Input file missing or empty: {merged_counts}")
            
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "BaseCellCalling", "BaseCellCalling.step1.py"),
            "--infile", merged_counts,
            "--outfile", output_prefix,
            "--ref", ref_genome
        ]
        print(f"Running variant calling step1 for {output_prefix}...")
        print(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Verify output
        if not os.path.exists(step1_output) or os.path.getsize(step1_output) == 0:
            raise RuntimeError(f"Output file missing or empty: {step1_output}")
            
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Variant calling step1 failed for {output_prefix}")
        print(f"Exit code: {e.returncode}")
        print(f"Error output:\n{e.stderr[:500]}...")
        return False
    except subprocess.TimeoutExpired:
        print(f"ERROR: Variant calling step1 timed out for {output_prefix}")
        return False
    except Exception as e:
        print(f"ERROR: {str(e)}")
        return False

def run_variant_calling_step2(args):
    """Run BaseCellCalling step 2 with enhanced validation"""
    step1_output, output_prefix, editing_sites, pon_file, scomatic_scripts_dir = args
    step2_output = f"{output_prefix}.calling.step2.tsv"
    
    try:
        # Verify step1 output
        if not os.path.exists(step1_output) or os.path.getsize(step1_output) == 0:
            raise ValueError(f"Step1 output missing or empty: {step1_output}")
            
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "BaseCellCalling", "BaseCellCalling.step2.py"),
            "--infile", step1_output,
            "--outfile", output_prefix,
            "--editing", editing_sites,
            "--pon", pon_file
        ]
        print(f"Running variant calling step2 for {output_prefix}...")
        print(f"Command: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd, 
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Verify output
        if not os.path.exists(step2_output) or os.path.getsize(step2_output) == 0:
            raise RuntimeError(f"Output file missing or empty: {step2_output}")
            
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Variant calling step2 failed for {output_prefix}")
        print(f"Exit code: {e.returncode}")
        print(f"Error output:\n{e.stderr[:500]}...")
        return False
    except subprocess.TimeoutExpired:
        print(f"ERROR: Variant calling step2 timed out for {output_prefix}")
        return False
    except Exception as e:
        print(f"ERROR: {str(e)}")
        return False

def filter_with_bed(input_tsv, bed_file, output_tsv):
    """Filter variants using bedtools intersect with validation"""
    try:
        # Check inputs
        if not os.path.exists(input_tsv) or os.path.getsize(input_tsv) == 0:
            print(f"ERROR: Input file missing or empty: {input_tsv}")
            return False
            
        if not os.path.exists(bed_file) or os.path.getsize(bed_file) == 0:
            print(f"ERROR: BED file missing or empty: {bed_file}")
            return False
            
        cmd = f"bedtools intersect -header -a {input_tsv} -b {bed_file} | " \
              f"awk '$1 ~ /^#/ || $6 == \"PASS\"' > {output_tsv}"
              
        print(f"Filtering {input_tsv} with BED file...")
        subprocess.run(cmd, shell=True, check=True)
        
        # Verify output
        if not os.path.exists(output_tsv) or os.path.getsize(output_tsv) == 0:
            print(f"WARNING: Output file empty after filtering: {output_tsv}")
            # Create placeholder
            with open(output_tsv, 'w') as f:
                f.write("## Empty after filtering\n")
                
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: BED filtering failed: {e}")
        return False
    except Exception as e:
        print(f"ERROR: Unexpected error in filter_with_bed: {str(e)}")
        return False

def prepare_sample_args(run_acc, adata_pp, working_dir, scomatic_scripts_dir, ref_genome):
    """Prepare arguments for sample processing with phase completion validation"""
    if not all([adata_pp, working_dir, scomatic_scripts_dir, ref_genome]):
        raise ValueError("Missing required arguments")
    
    # Get sample metadata from AnnData
    mask = adata_pp.obs['run_accession'] == run_acc
    if not any(mask):
        print(f"SKIPPING {run_acc}: Not found in AnnData")
        return None
        
    subset = adata_pp[mask, :]
    series_id = subset.obs['series_id'].iloc[0]
    
    # Define directory paths
    scomatic_dir = os.path.join(working_dir, 'fastq', series_id, run_acc, f"{run_acc}_S1_L001_", "outs", "SComatic")
    variant_calling_dir = os.path.join(scomatic_dir, 'VariantCalling')
    
    # ===== PHASE COMPLETION VALIDATION =====
    
    # 1. Check if Phase 1 completed (Annotation file exists)
    meta_file = os.path.join(scomatic_dir, 'cell_barcode_annotation.tsv')
    if not os.path.exists(meta_file):
        print(f"SKIPPING {run_acc}: Missing annotation file (Phase 1 failed)")
        return None
    
    # 2. Check if Phase 1 completed (Split BAM files exist)
    bam_files = glob.glob(os.path.join(scomatic_dir, f"{run_acc}.*.bam"))
    if not bam_files:
        print(f"SKIPPING {run_acc}: No split BAM files found (Phase 1 failed)")
        return None
    
    # 3. Check if Phase 2 completed (BaseCellCounts directory exists with files)
    basecell_counts_dir = os.path.join(scomatic_dir, 'BaseCellCounts')
    if not os.path.exists(basecell_counts_dir):
        print(f"SKIPPING {run_acc}: BaseCellCounts directory missing (Phase 2 failed)")
        return None
    
    basecell_files = glob.glob(os.path.join(basecell_counts_dir, "*.tsv"))
    if not basecell_files:
        print(f"SKIPPING {run_acc}: No BaseCellCount files found (Phase 2 failed)")
        return None
    
    # 4. Check if Phase 3 completed (Merged counts file exists)
    merged_counts_file = os.path.join(scomatic_dir, 'MergedCounts', f"{run_acc}.BaseCellCounts.AllCellTypes.tsv")
    if not os.path.exists(merged_counts_file) or os.path.getsize(merged_counts_file) == 0:
        print(f"SKIPPING {run_acc}: Merged counts file missing or empty (Phase 3 failed)")
        return None
    
    # 5. Check if Phase 4 completed (Variant calling step2 file exists)
    step2_file = os.path.join(variant_calling_dir, f"{run_acc}.calling.step2.tsv")
    if not os.path.exists(step2_file) or os.path.getsize(step2_file) == 0:
        print(f"SKIPPING {run_acc}: Variant calling step2 file missing (Phase 4 failed)")
        return None
    
    # 6. Check if Phase 5 completed (Filtered variants exist)
    filtered_dir = os.path.join(scomatic_dir, 'FilteredVariants')
    filtered_files = glob.glob(os.path.join(filtered_dir, "*.filtered.tsv"))
    if not filtered_files:
        print(f"SKIPPING {run_acc}: No filtered variant files found (Phase 5 failed)")
        return None
    
    # ===== PREPARE PHASE 6 ARGUMENTS =====
    
    # 1. Callable sites per cell type
    callable_args = None
    callable_sites_dir = os.path.join(scomatic_dir, 'CellTypeCallableSites')
    step1_file = os.path.join(variant_calling_dir, f"{run_acc}.calling.step1.tsv")
    
    if os.path.exists(step1_file) and os.path.getsize(step1_file) > 0:
        os.makedirs(callable_sites_dir, exist_ok=True)
        callable_args = (step1_file, os.path.join(callable_sites_dir, run_acc), "150", "2", scomatic_scripts_dir)
    else:
        print(f"WARNING: Step1 file missing for callable sites: {step1_file}")
    
    # 2. Sites per cell
    sites_per_cell_args = []
    sites_per_cell_dir = os.path.join(scomatic_dir, 'UniqueCellCallableSites')
    
    if os.path.exists(step1_file) and os.path.getsize(step1_file) > 0:
        for bam in bam_files:
            if os.path.exists(bam) and os.path.getsize(bam) > 0:
                temp_dir = os.path.join(sites_per_cell_dir, f"temp_{os.path.basename(bam).split('.')[-2]}")
                sites_per_cell_args.append((bam, step1_file, ref_genome, sites_per_cell_dir, temp_dir, scomatic_scripts_dir))
    
    # 3. Single cell genotypes
    genotype_args = []
    single_cell_dir = os.path.join(scomatic_dir, 'SingleCellAlleles')
    
    if os.path.exists(step2_file) and os.path.getsize(step2_file) > 0:
        for bam in bam_files:
            cell_type = os.path.basename(bam).split('.')[-2]
            output_file = os.path.join(single_cell_dir, f"{cell_type}.single_cell_genotype.tsv")
            temp_dir = os.path.join(single_cell_dir, f"temp_{cell_type}")
            genotype_args.append((bam, step2_file, meta_file, ref_genome, output_file, temp_dir, scomatic_scripts_dir))
    
    # 4. Filtered genotypes
    filtered_genotype_files = []
    single_cell_dir = os.path.join(scomatic_dir, 'SingleCellAlleles')
    if os.path.exists(step2_file) and os.path.getsize(step2_file) > 0:
        filtered_genotype_files = [(step2_file, single_cell_dir)]
    
    # ===== RETURN ARGUMENTS =====
    
    return {
        'run_acc': run_acc,
        'callable_args': callable_args,
        'sites_per_cell_args': sites_per_cell_args,
        'genotype_args': genotype_args,
        'filtered_genotype_files': filtered_genotype_files,
        'step1_file': step1_file,
        'directories': {
            'root': scomatic_dir,
            'callable_sites': callable_sites_dir,
            'sites_per_cell': sites_per_cell_dir,
            'single_cell': single_cell_dir,
            'variant_calling': variant_calling_dir
        }
    }   

def run_callable_sites(args):
    """Compute callable sites per cell type"""
    step1_file, output_prefix, max_cov, min_cell_types, scomatic_scripts_dir = args
    try:
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "GetCallableSites", "GetAllCallableSites.py"),
            "--infile", step1_file,
            "--outfile", output_prefix,
            "--max_cov", str(max_cov),
            "--min_cell_types", str(min_cell_types)
        ]
        print(f"Computing callable sites for {output_prefix}...")
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Callable sites computation failed: {e}")
        return False

def run_sites_per_cell(args):
    """Compute callable sites per cell with robust error handling"""
    bam_file, step1_file, ref_genome, output_dir, temp_dir, scomatic_scripts_dir = args
    cell_type = os.path.basename(bam_file).split('.')[-2]
    output_file = os.path.join(output_dir, f"{os.path.basename(bam_file)}.SitesPerCell.tsv")
    
    try:
        # Validate inputs
        if not os.path.exists(bam_file) or os.path.getsize(bam_file) == 0:
            print(f"WARNING: BAM file missing or empty: {bam_file}")
            return False
            
        if not os.path.exists(step1_file) or os.path.getsize(step1_file) == 0:
            print(f"WARNING: Step1 file missing or empty: {step1_file}")
            return False

        # Create directories
        os.makedirs(temp_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)

        # Remove existing output file if present
        if os.path.exists(output_file):
            print(f"Removing existing output file: {output_file}")
            try:
                os.remove(output_file)
            except Exception as e:
                print(f"WARNING: Could not remove existing file - {e}")
        
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "SitesPerCell", "SitesPerCell.py"),
            "--bam", bam_file,
            "--infile", step1_file,
            "--ref", ref_genome,
            "--out_folder", output_dir,
            "--tmp_dir", temp_dir,
            "--nprocs", "1"
        ]
        print(f"\nProcessing sites per cell for {cell_type}...")
        print(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd, 
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
        except Exception as e:
            print(f"ERROR: An error occured {e}")
            return False
        
        # Check output
        if not os.path.exists(output_file):
            print(f"ERROR: Output file not created: {output_file}")
            return False
            
        # Handle empty output
        if os.path.getsize(output_file) < 20:  # Less than header size
            print(f"WARNING: Output file empty, creating placeholder: {output_file}")
            with open(output_file, 'w') as f:
                f.write("CB,SitesPerCell\n")
                
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command failed for {cell_type}")
        print(f"Exit code: {e.returncode}")
        print(f"Error output:\n{e.stderr[:500]}...")
        return False
    except Exception as e:
        print(f"ERROR: Processing failed for {cell_type}: {str(e)}")
        return False
    finally:
        # Keep temp dirs for debugging failures
        if 'result' in locals() and result.returncode != 0 and os.path.exists(temp_dir):
            print(f"Keeping temp dir for debugging: {temp_dir}")

def run_single_cell_genotype(args):
    """Compute genotypes for each cell with enhanced validation"""
    bam_file, variant_file, meta_file, ref_genome, output_file, temp_dir, scomatic_scripts_dir = args
    cell_type = os.path.basename(bam_file).split('.')[-2]
    
    try:
        # Add debug logging
        print(f"Processing genotype for: {os.path.basename(bam_file)}")
        print(f"Using meta file: {meta_file}")
        print(f"Output will be: {output_file}")
        
        # Check if output already exists
        if os.path.exists(output_file):
            print(f"Output file already exists, skipping: {output_file}")
            return True
            
        # Validate inputs
        if not os.path.exists(bam_file) or os.path.getsize(bam_file) == 0:
            print(f"ERROR: BAM file missing or empty: {bam_file}")
            return False
            
        if not os.path.exists(variant_file) or os.path.getsize(variant_file) == 0:
            print(f"ERROR: Variant file missing or empty: {variant_file}")
            return False
            
        # CRITICAL: Validate meta file exists
        if not os.path.exists(meta_file) or os.path.getsize(meta_file) == 0:
            print(f"ERROR: Meta file missing or empty: {meta_file}")
            return False
            
        os.makedirs(temp_dir, exist_ok=True)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "SingleCellGenotype", "SingleCellGenotype.py"),
            "--bam", bam_file,
            "--infile", variant_file,
            "--nprocs", "1",
            "--meta", meta_file,
            "--outfile", output_file,
            "--tmp_dir", temp_dir,
            "--ref", ref_genome
        ]
        print(f"Computing genotypes for {cell_type}...")
        print(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd, 
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            print(f"Successfully processed genotype for {cell_type}")
        except Exception as e:
            print(f"ERROR: An error occurred {e}")
            if hasattr(e, 'stderr') and e.stderr:
                print(f"Error output: {e.stderr[:500]}...")
            return False
        
        # Verify output
        if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            print(f"ERROR: Output file not created or empty: {output_file}")
            return False
            
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Single cell genotype failed for {bam_file}: {e}")
        print(f"Error output:\n{e.stderr[:500]}...")
        return False
    finally:
        # Clean up temp dir only on success
        if os.path.exists(temp_dir) and 'result' in locals() and result.returncode == 0:
            shutil.rmtree(temp_dir, ignore_errors=True)

def filter_and_annotate_sc_genotypes(input_dir, variant_file, output_root, sample_id):
    """Filter single cell genotypes and add trinucleotide context"""
    # Create output directory structure
    final_output_dir = os.path.join(output_root, 'FilteredSingleCellAlleles')
    os.makedirs(final_output_dir, exist_ok=True)
    
    # Load variant information with corrected base extraction
    print(f"Loading variant context information from {variant_file}...")
    variants_dict = {}
    context_issues = 0
    
    try:
        with open(variant_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or line.startswith('Chr'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    context_issues += 1
                    continue
                    
                chrom = parts[0]
                pos = parts[1]
                ref = parts[3]
                alt = parts[4]
                up_context = parts[7]  # Up_context column
                down_context = parts[8]  # Down_context column
                
                # OLD REMOVED 10-02-2025 
                # Debug print for first 5 variants
                # if line_num <= 5:
                #     print(f"Variant {line_num}: {chrom}:{pos} {ref}>{alt}")
                #     print(f"  Up_context: '{up_context}'")
                #     print(f"  Down_context: '{down_context}'")
                
                # CORRECTED BASE EXTRACTION:
                # Upstream base is the SECOND character (index 1) of Up_context
                # Downstream base is the SECOND character (index 1) of Down_context
                up_base = up_context[1] if len(up_context) >= 2 else 'N'
                down_base = down_context[1] if len(down_context) >= 2 else 'N'
                
                # if line_num <= 5:
                #     print(f"  Extracted bases: upstream='{up_base}' downstream='{down_base}'")
                #     print(f"  Expected trinucleotide: {up_base}{ref}{down_base} -> {up_base}{alt}{down_base}")
                
                key = f"{chrom}_{pos}"
                variants_dict[key] = {
                    'ref': ref,
                    'alt': alt,
                    'up_base': up_base,
                    'down_base': down_base
                }
                
    except Exception as e:
        print(f"Error loading variant file {variant_file}: {str(e)}")
        return [], None
    
    if context_issues > 0:
        print(f"Warning: {context_issues} variants had incomplete context information")

    # Process each cell type file
    print(f"Processing cell type files in {input_dir}...")
    input_files = glob.glob(os.path.join(input_dir, '*.single_cell_genotype.tsv'))
    
    if not input_files:
        print("Warning: No input files found matching pattern *.single_cell_genotype.tsv")
        return [], None
    
    processed_files = []
    all_cells_data = []

    # Tracking statistics
    total_rows = 0
    filtered_by_base_mismatch = 0
    filtered_by_alt_reads = 0
    filtered_by_total_depth = 0
    rows_kept = 0
    
    for input_file in input_files:
        try:
            # Skip empty files
            if os.path.getsize(input_file) == 0:
                print(f"WARNING: Skipping empty file: {input_file}")
                continue
                
            df = pd.read_csv(input_file, sep='\t')
            if df.empty:
                print(f"WARNING: Empty dataframe in {input_file}")
                continue
            
            initial_count = len(df)
            total_rows += initial_count
            
            def get_alt_bases(alt_str):
                """Extract unique bases from ALT_expected, handling commas and pipes"""
                if pd.isna(alt_str) or alt_str == '.':
                    return set()
                # Replace pipes with commas, split, and get unique bases
                alt_str = str(alt_str).replace('|', ',')
                bases = [base.strip() for base in alt_str.split(',') if base.strip()]
                return set(bases)
            
            df['ALT_bases_set'] = df['ALT_expected'].apply(get_alt_bases)
            
            # Filter: Base_observed must be in the set of ALT bases (not just first character)
            def base_matches_alt(row):
                """Check if Base_observed matches any of the ALT_expected bases"""
                if not row['ALT_bases_set']:
                    return False
                return row['Base_observed'] in row['ALT_bases_set']
            
            mask_base = df.apply(base_matches_alt, axis=1)
            filtered_by_base_mismatch += (~mask_base).sum()
            
            # Filter 2: Must have >= 3 reads supporting ALT
            mask_alt_reads = df['Num_reads'] >= 3
            filtered_by_alt_reads += (~mask_alt_reads & mask_base).sum()
            
            # Filter 3: Must have >= 5 total depth
            mask_total_depth = df['Total_depth'] >= 5
            filtered_by_total_depth += (~mask_total_depth & mask_base & mask_alt_reads).sum()
            
            # Combined filter
            df_filtered = df[mask_base & mask_alt_reads & mask_total_depth].copy()
            rows_kept += len(df_filtered)
            
            # Drop the temporary column
            df_filtered = df_filtered.drop(columns=['ALT_bases_set'])
            
            if df_filtered.empty:
                print(f"No matching variants after filtering in {input_file}")
                continue
                
            # Add barcode suffix
            df_filtered['CB'] = df_filtered['CB'] + f"-1-{sample_id}"
            
            # Add trinucleotide context
            def get_trinucleotide(row):
                chrom_pos = f"{row['#CHROM']}_{row['Start']}"
                variant_info = variants_dict.get(chrom_pos, {})
                ref_tri = variant_info.get('up_base', 'N') + row['REF'] + variant_info.get('down_base', 'N')
                alt_tri = variant_info.get('up_base', 'N') + row['Base_observed'] + variant_info.get('down_base', 'N')
                return pd.Series([ref_tri, alt_tri])
            
            df_filtered[['REF_TRI', 'ALT_TRI']] = df_filtered.apply(get_trinucleotide, axis=1)
            
            # Save cell type file
            cell_type = os.path.basename(input_file).split('.')[0]
            output_file = os.path.join(final_output_dir, f"{cell_type}.single_cell_genotype.filtered.tsv")
            df_filtered.to_csv(output_file, sep='\t', index=False)
            processed_files.append(output_file)
            
            # Add to combined data
            all_cells_data.append(df_filtered)
            
        except Exception as e:
            print(f"ERROR processing {input_file}: {str(e)}")
            continue
    
    # Print filtering statistics
    print("\n" + "="*60)
    print("FILTERING STATISTICS")
    print("="*60)
    print(f"Total rows processed: {total_rows}")
    print(f"Filtered by base mismatch (Base_observed != ALT_expected): {filtered_by_base_mismatch}")
    print(f"Filtered by insufficient ALT reads (<3): {filtered_by_alt_reads}")
    print(f"Filtered by insufficient total depth (<5X): {filtered_by_total_depth}")
    print(f"Rows kept after all filters: {rows_kept}")
    print(f"Percentage retained: {100 * rows_kept / total_rows:.2f}%")
    print("="*60 + "\n")
    
    # Create combined output
    if all_cells_data:
        combined_cells = pd.concat(all_cells_data, ignore_index=True)
        combined_output = os.path.join(final_output_dir, "all_cell.single_cell_genotype.filtered.tsv")
        combined_cells.to_csv(combined_output, sep='\t', index=False)
        processed_files.append(combined_output)
        return processed_files, combined_cells
    return processed_files, None

def run_trinucleotide_context(args):
    """Compute trinucleotide context background"""
    step1_files_list, output_file, scomatic_scripts_dir = args
    try:
        # Create temporary directory if it doesn't exist
        temp_dir = os.path.join(os.path.dirname(output_file), "temp_trinuc")
        os.makedirs(temp_dir, exist_ok=True)
        
        # Create a temporary input file list
        input_list_file = os.path.join(temp_dir, "trinuc_input_list.txt")
        
        # Write all input files to the temporary list file
        with open(input_list_file, 'w') as f:
            for file in step1_files_list:
                if os.path.exists(file):  # Only include files that exist
                    f.write(f"{file}\n")
        
        if os.path.getsize(input_list_file) == 0:
            print("WARNING: No valid input files for trinucleotide context")
            return False
        
        cmd = [
            "python",
            os.path.join(scomatic_scripts_dir, "TrinucleotideBackground", "TrinucleotideContextBackground.py"),
            "--in_tsv", input_list_file,
            "--out_file", output_file
        ]
        print("Computing trinucleotide context background...")
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Trinucleotide context failed: {e}")
        return False
    except Exception as e:
        print(f"ERROR: Unexpected error in trinucleotide context: {str(e)}")
        return False
    finally:
        # Clean up the temporary directory
        if 'temp_dir' in locals() and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir, ignore_errors=True)

def generate_complete_callable_sites(working_dir, valid_samples, adata_pp, cell_annotations):
    """Generate complete callable sites file including all cells with diagnostics"""
    combined_callable_dir = os.path.join(working_dir, 'CombinedCallableSites')
    os.makedirs(combined_callable_dir, exist_ok=True)
    
    # Create complete cell list
    all_cells = pd.DataFrame({'CB': cell_annotations['cell_barcodes']})
    all_cells['SitesPerCell'] = 0
    
    # Track processing stats
    stats = {
        'total_cells': len(all_cells),
        'processed_cells': 0,
        'missing_cells': 0,
        'samples_processed': 0,
        'samples_missing': 0,
        'files_processed': 0,
        'files_skipped': 0
    }
    
    for run_acc in tqdm(valid_samples, desc="Processing callable sites"):
        mask = adata_pp.obs['run_accession'] == run_acc
        subset = adata_pp[mask, :]
        series_id = subset.obs['series_id'].iloc[0]
        
        callable_sites_dir = os.path.join(
            working_dir, 'fastq', series_id, run_acc,
            f"{run_acc}_S1_L001_", "outs", "SComatic", 'UniqueCellCallableSites'
        )
        
        if not os.path.exists(callable_sites_dir):
            stats['samples_missing'] += 1
            print(f"Sample {run_acc}: Callable sites directory missing: {callable_sites_dir}")
            continue
            
        stats['samples_processed'] += 1
        
        # Look for the specific file pattern generated by run_sites_per_cell
        callable_files = glob.glob(os.path.join(callable_sites_dir, "*.SitesPerCell.tsv"))
        
        if not callable_files:
            # Fallback to any TSV files if the specific pattern isn't found
            callable_files = glob.glob(os.path.join(callable_sites_dir, "*.tsv"))
            print(f"Sample {run_acc}: Using fallback file pattern, found {len(callable_files)} files")
        
        for file_path in callable_files:
            try:
                # Skip small files (less than header size)
                if os.path.getsize(file_path) < 20:
                    stats['files_skipped'] += 1
                    continue
                
                # Try to detect the delimiter automatically
                try:
                    with open(file_path, 'r') as f:
                        first_line = f.readline()
                    
                    if ',' in first_line:
                        df = pd.read_csv(file_path, sep=',')
                    elif '\t' in first_line:
                        df = pd.read_csv(file_path, sep='\t')
                    else:
                        print(f"WARNING: Could not detect delimiter in {file_path}, skipping")
                        stats['files_skipped'] += 1
                        continue
                        
                except Exception as e:
                    print(f"Error reading {file_path}: {str(e)}")
                    stats['files_skipped'] += 1
                    continue
                
                if df.empty:
                    stats['files_skipped'] += 1
                    continue
                
                # Check if the file has the expected columns
                if 'CB' not in df.columns or 'SitesPerCell' not in df.columns:
                    print(f"WARNING: File {os.path.basename(file_path)} missing required columns, skipping")
                    stats['files_skipped'] += 1
                    continue
                    
                # Format CB with sample identifier
                df['CB'] = df['CB'] + f"-1-{run_acc}"
                
                # Update sites count for matching cells
                cells_updated = 0
                for _, row in df.iterrows():
                    if row['CB'] in all_cells['CB'].values:
                        idx = all_cells.index[all_cells['CB'] == row['CB']][0]
                        all_cells.at[idx, 'SitesPerCell'] = row['SitesPerCell']
                        cells_updated += 1
                
                stats['processed_cells'] += cells_updated
                stats['files_processed'] += 1
                
                if cells_updated == 0:
                    print(f"WARNING: No cells matched from {os.path.basename(file_path)}")
                    
            except Exception as e:
                print(f"Error processing {file_path}: {str(e)}")
                stats['files_skipped'] += 1
                continue
            
    # Calculate missing cells
    stats['missing_cells'] = stats['total_cells'] - stats['processed_cells']
    
    # Save output
    output_path = os.path.join(combined_callable_dir, 'complete_callable_sites.tsv')
    all_cells.to_csv(output_path, sep='\t', index=False)
    
    # Save diagnostics
    stats_path = os.path.join(combined_callable_dir, 'callable_sites_stats.txt')
    with open(stats_path, 'w') as f:
        f.write("Callable Sites Processing Statistics\n")
        f.write("="*50 + "\n")
        f.write(f"Total cells: {stats['total_cells']}\n")
        f.write(f"Processed cells: {stats['processed_cells']}\n")
        f.write(f"Missing cells: {stats['missing_cells']}\n")
        f.write(f"Samples processed: {stats['samples_processed']}\n")
        f.write(f"Samples missing: {stats['samples_missing']}\n")
        f.write(f"Files processed: {stats['files_processed']}\n")
        f.write(f"Files skipped: {stats['files_skipped']}\n")
        f.write(f"Success rate: {stats['files_processed']/(stats['files_processed'] + stats['files_skipped'])*100:.1f}%\n")
    
    print(f"Saved complete callable sites to: {output_path}")
    print(f"Statistics: {stats}")
    
    # Additional debug output
    print("Callable sites summary:")
    print(f"  Total cells in annotation: {stats['total_cells']}")
    print(f"  Cells with callable sites data: {stats['processed_cells']}")
    print(f"  Cells missing data: {stats['missing_cells']}")
    print(f"  Samples processed: {stats['samples_processed']}")
    print(f"  Samples missing: {stats['samples_missing']}")
    
    return output_path

def process_sample(sample_args):
    """Process all steps for a single sample"""
    from os import cpu_count
    run_acc = sample_args['run_acc']
    results = {'run_acc': run_acc, 'success': True}
    
    try:
        # Create directories
        for dir_path in sample_args['directories'].values():
            if dir_path:
                os.makedirs(dir_path, exist_ok=True)
        
        # 1. Callable sites per cell type
        if sample_args['callable_args']:
            try:
                print(f"Processing callable sites for {run_acc}...")
                success = run_callable_sites(sample_args['callable_args'])
                results['callable_sites'] = success
                if not success:
                    results['success'] = False
                    print(f"Callable sites failed for {run_acc}")
            except Exception as e:
                print(f"ERROR in callable sites for {run_acc}: {str(e)}")
                results['callable_sites'] = False
                results['success'] = False
        
        # 2. Sites per cell
        if sample_args['sites_per_cell_args']:
            try:
                print(f"Processing sites per cell for {run_acc}...")
                with ThreadPool(min(cpu_count(), len(sample_args['sites_per_cell_args']))) as pool:
                    sites_results = list(pool.imap(run_sites_per_cell, sample_args['sites_per_cell_args']))
                results['sites_per_cell'] = sum(sites_results)
                if sum(sites_results) < len(sites_results):
                    results['success'] = False
                    print(f"Sites per cell partially failed for {run_acc}: {sum(sites_results)}/{len(sites_results)} successful")
            except Exception as e:
                print(f"ERROR in sites per cell for {run_acc}: {str(e)}")
                results['sites_per_cell'] = 0
                results['success'] = False
        
        # 3. Single cell genotypes - CRITICAL SECTION WITH FIXES
        if sample_args['genotype_args']:
            try:
                print(f"Processing single cell genotypes for {run_acc}...")
                print(f"Number of genotype tasks: {len(sample_args['genotype_args'])}")
                
                # Debug: print first few genotype args to verify inputs
                for i, args in enumerate(sample_args['genotype_args'][:3]):
                    print(f"Genotype task {i}: BAM={os.path.basename(args[0])}, Meta={os.path.basename(args[2])}")
                
                with ThreadPool(min(cpu_count(), len(sample_args['genotype_args']))) as pool:
                    genotype_results = list(pool.imap(run_single_cell_genotype, sample_args['genotype_args']))
                
                results['single_cell_genotypes'] = sum(genotype_results)
                
                # CRITICAL FIX: If any genotype fails, skip filtered genotypes step
                if sum(genotype_results) < len(genotype_results):
                    results['success'] = False
                    failed_count = len(genotype_results) - sum(genotype_results)
                    print(f"WARNING: {failed_count} genotype processes failed for {run_acc}")
                    # Skip filtered genotypes since inputs are incomplete
                    sample_args['filtered_genotype_files'] = []
                    results['filtered_genotypes'] = False
                    
            except Exception as e:
                print(f"ERROR in genotype processing for {run_acc}: {str(e)}")
                import traceback
                traceback.print_exc()
                results['single_cell_genotypes'] = 0
                results['success'] = False
                sample_args['filtered_genotype_files'] = []  # Skip next step
                results['filtered_genotypes'] = False
        
        # 4. Filtered genotypes - Only run if previous steps were successful
        if sample_args['filtered_genotype_files'] and results.get('single_cell_genotypes', 0) > 0:
            try:
                print(f"Processing filtered genotypes for {run_acc}...")
                variant_file, input_dir = sample_args['filtered_genotype_files'][0]
                
                # Verify inputs exist before proceeding
                if not os.path.exists(variant_file) or os.path.getsize(variant_file) == 0:
                    print(f"ERROR: Variant file missing or empty: {variant_file}")
                    raise FileNotFoundError(f"Variant file not found: {variant_file}")
                
                if not os.path.exists(input_dir):
                    print(f"ERROR: Input directory missing: {input_dir}")
                    raise FileNotFoundError(f"Input directory not found: {input_dir}")
                
                processed_files, combined_data = filter_and_annotate_sc_genotypes(
                    input_dir, variant_file, sample_args['directories']['root'], run_acc
                )
                results['filtered_genotypes'] = bool(processed_files)
                results['combined_data'] = combined_data
                
                if processed_files:
                    print(f"Successfully processed {len(processed_files)} filtered genotype files for {run_acc}")
                else:
                    print(f"No filtered genotype files produced for {run_acc}")
                    
            except Exception as e:
                print(f"ERROR in genotype filtering for {run_acc}: {str(e)}")
                import traceback
                traceback.print_exc()
                results['filtered_genotypes'] = False
                results['success'] = False
        else:
            # Log why filtered genotypes were skipped
            if not sample_args['filtered_genotype_files']:
                print(f"Filtered genotypes skipped for {run_acc}: No genotype files specified")
            elif results.get('single_cell_genotypes', 0) == 0:
                print(f"Filtered genotypes skipped for {run_acc}: No successful genotype processing")
            results['filtered_genotypes'] = False
                
        return results
        
    except Exception as e:
        print(f"CRITICAL ERROR processing sample {run_acc}: {str(e)}")
        import traceback
        traceback.print_exc()
        return {'run_acc': run_acc, 'success': False, 'error': str(e)}

def monitor_resources(interval=300):
    """Monitor system resources during processing"""
    while not monitor_resources.stop:
        cpu = psutil.cpu_percent()
        mem = psutil.virtual_memory().percent
        disk = psutil.disk_usage('/').percent
        print(f"Resource Monitor | CPU: {cpu}% | Mem: {mem}% | Disk: {disk}%")
        time.sleep(interval)

def verify_pipeline_output(final_mutations_path, callable_sites_path, cell_annotations_path, trinuc_path=None):
    """Validate final pipeline output integrity"""
    print("\n" + "="*80)
    print("FINAL OUTPUT VERIFICATION")
    print("="*80)
    
    # Load data
    try:
        cell_annotations = pd.read_csv(cell_annotations_path, sep='\t')
        expected_cells = set(cell_annotations['cell_barcodes'])
    except Exception as e:
        print(f"ERROR loading cell annotations: {str(e)}")
        return
        
    try:
        if os.path.exists(final_mutations_path):
            mut_df = pd.read_csv(final_mutations_path, sep='\t')
            mutation_cells = set(mut_df['CB']) if 'CB' in mut_df.columns else set()
        else:
            mutation_cells = set()
            print("ERROR: Mutations file not found")
    except Exception as e:
        print(f"ERROR loading mutations: {str(e)}")
        mutation_cells = set()
        
    try:
        if os.path.exists(callable_sites_path):
            callable_df = pd.read_csv(callable_sites_path, sep='\t')
            callable_cells = set(callable_df['CB']) if 'CB' in callable_df.columns else set()
        else:
            callable_cells = set()
            print("ERROR: Callable sites file not found")
    except Exception as e:
        print(f"ERROR loading callable sites: {str(e)}")
        callable_cells = set()
    
    if trinuc_path:
        try:
            if os.path.exists(trinuc_path):
                trinuc_df = pd.read_csv(trinuc_path, sep='\t')
                print(f"Trinucleotide background file valid: {trinuc_path}")
                print(f"Contains {len(trinuc_df)} entries")
            else:
                print(f"WARNING: Trinucleotide file not found: {trinuc_path}")
        except Exception as e:
            print(f"ERROR verifying trinucleotide file: {str(e)}")    
            
    # Generate report
    report = [
        "Pipeline Output Verification Report",
        "="*50,
        f"Expected cells: {len(expected_cells)}",
        f"Cells in mutations file: {len(mutation_cells)}",
        f"Cells in callable sites: {len(callable_cells)}",
        "",
        f"Cells missing from mutations: {len(expected_cells - mutation_cells)}",
        f"Cells missing from callable sites: {len(expected_cells - callable_cells)}",
        f"Cells with mutations but no callable sites: {len(mutation_cells - callable_cells)}"
    ]
    
    # Print and save report
    print("\n".join(report))
    report_path = os.path.join(os.path.dirname(final_mutations_path), "pipeline_verification_report.txt")
    with open(report_path, 'w') as f:
        f.write("\n".join(report))
    print(f"Verification report saved to: {report_path}")
    


#%%
#=============================#
#        Main Pipeline        #
#=============================#
def main():
    import os
    import yaml
    import pandas as pd
    import scanpy as sc
    from tqdm import tqdm
    import pysam
    import subprocess
    import shutil
    import glob
    from multiprocessing.pool import ThreadPool
    from math import ceil
    from time import perf_counter
    from os import cpu_count
    import pickle

    # Start resource monitoring
    monitor_resources.stop = False
    monitor_thread = threading.Thread(target=monitor_resources)
    monitor_thread.daemon = True
    monitor_thread.start()
    
    #TMP section to researt from a certain step
    config_path = '/work/sdz852/WORKING/SC/fastq/Breast_cancer/config_SABCS.yaml'
    with open(config_path, 'r') as config:
        config_yaml = yaml.safe_load(config)
    working_dir = config_yaml['output_dir']
    os.chdir(working_dir)
 
    RESTART_FROM_STEP = 1  # Set to 5 to restart from BED filtering
    CHECKPOINT_PATH = os.path.join(working_dir, "basecell_counter_args.pkl")
    
    # Initialize state variables we'll need for restart
    final_valid_samples = None
    adata_pp = None
    all_sample_args = []
    
    if RESTART_FROM_STEP >= 6:
        print(f"\n{'='*80}")
        print(f"RESTARTING PIPELINE FROM STEP {RESTART_FROM_STEP}")
        print(f"{'='*80}")
        
        if os.path.exists(CHECKPOINT_PATH):
            try:
                # Load Step 2 checkpoint
                with open(CHECKPOINT_PATH, 'rb') as f:
                    basecell_counter_args, basecell_results = pickle.load(f)
                
                # Reconstruct sample list from BAM paths
                run_accessions_step2 = set()
                for (bam_path, _, _, _) in basecell_counter_args:
                    # Extract run accession from BAM path
                    # Example: /path/to/SRR12345678.T_cell.bam
                    filename = os.path.basename(bam_path)
                    run_acc = filename.split('.')[0]
                    run_accessions_step2.add(run_acc)
                
                # Load AnnData object
                print("Reloading AnnData object...")
                adata_pp = sc.read_h5ad(os.path.join(working_dir, 'adata_pp.h5ad'))
                
                # Verify which samples completed Step 4
                final_valid_samples = []
                for run_acc in tqdm(run_accessions_step2, desc="Validating Step 4 completion"):
                    # Get sample metadata from AnnData
                    mask = adata_pp.obs['run_accession'] == run_acc
                    if not any(mask):
                        continue
                    subset = adata_pp[mask, :]
                    series_id = subset.obs['series_id'].iloc[0]
                    
                    # Check if Phase 5 completed (both filtered files AND step2.pass file exist)
                    variant_calling_dir = os.path.join(
                        working_dir, 'fastq', series_id, run_acc,
                        f"{run_acc}_S1_L001_", "outs", "SComatic", "VariantCalling"
                    )
                    filtered_dir = os.path.join(
                        working_dir, 'fastq', series_id, run_acc,
                        f"{run_acc}_S1_L001_", "outs", "SComatic", "FilteredVariants"
                    )
                    
                    # Check for both the filtered files AND the step2.pass file that Phase 6 needs
                    step2_file = os.path.join(variant_calling_dir, f"{run_acc}.calling.step2.tsv")
                    filtered_files = glob.glob(os.path.join(filtered_dir, "*.filtered.tsv"))
                    
                    if filtered_files and os.path.exists(step2_file) and os.path.getsize(step2_file) > 0:
                        final_valid_samples.append(run_acc)
                        print(f"Sample {run_acc} completed Phase 5 successfully")
                    else:
                        print(f"Sample {run_acc} skipped: filtered_files={bool(filtered_files)}, step2_pass_exists={os.path.exists(step2_file)}")
                
                print(f"Resuming with {len(final_valid_samples)} samples that completed Step 5")
                
                # Load configuration
                config_path = '/work/sdz852/WORKING/SC/fastq/Breast_cancer/config_SABCS.yaml'
                with open(config_path, 'r') as config:
                    config_yaml = yaml.safe_load(config)
                working_dir = config_yaml['output_dir']
                os.chdir(working_dir)
                
                # Define paths
                scomatic_scripts_dir = '/work/sdz852/WORKING/SComatic/scripts'
                ref_genome = '/work/sdz852/WORKING/SC/ref/GRCh38/fasta/genome.fa'
                editing_sites = '/work/sdz852/WORKING/SComatic/RNAediting/AllEditingSites.hg38.txt'
                pon_file = '/work/sdz852/WORKING/SComatic/PoNs/PoN.scRNAseq.hg38.tsv'
                bed_file = '/work/sdz852/WORKING/SComatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed'    
                NCPUS = cpu_count()  # Limit to 16 cores to avoid overloading
                
                # Load processed AnnData object
                print("Loading AnnData object...")
                adata_pp = sc.read_h5ad(os.path.join(working_dir, 'adata_pp.h5ad'))
                
                # Rebuild sample arguments for Step 6
                if RESTART_FROM_STEP == 6:
                    print("Preparing sample arguments for Step 6...")
                    all_sample_args = [prepare_sample_args(run_acc, adata_pp, working_dir, scomatic_scripts_dir, ref_genome) for run_acc in final_valid_samples]
                
            except Exception as e:
                print(f"Error loading Step 2 checkpoint: {str(e)}")
                import traceback
                traceback.print_exc()
                return
        else:
            print(f"Step 2 checkpoint not found: {CHECKPOINT_PATH}")
            return     
    
    
    try:
        # Load configuration
        config_path = '/work/sdz852/WORKING/SC/fastq/Breast_cancer/config_SABCS.yaml'
        with open(config_path, 'r') as config:
            config_yaml = yaml.safe_load(config)
        working_dir = config_yaml['output_dir']
        os.chdir(working_dir)
        
        # Define paths
        scomatic_scripts_dir = '/work/sdz852/WORKING/SComatic/scripts'
        ref_genome = '/work/sdz852/WORKING/SC/ref/GRCh38/fasta/genome.fa'
        editing_sites = '/work/sdz852/WORKING/SComatic/RNAediting/AllEditingSites.hg38.txt'
        pon_file = '/work/sdz852/WORKING/SComatic/PoNs/PoN.scRNAseq.hg38.tsv'
        bed_file = '/work/sdz852/WORKING/SComatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed'
        # For phase 8
        final_output_path = os.path.join(working_dir, 'all_samples.single_cell_genotype.filtered.tsv')
        callable_path = os.path.join(working_dir, 'CombinedCallableSites', 'complete_callable_sites.tsv')
        cell_annotations_path = os.path.join(working_dir, 'cell_annotations.txt')
        output_trinuc = os.path.join(working_dir, 'trinucleotide_background.tsv')
        NCPUS = cpu_count()  # Limit to 16 cores to avoid overloading
        
        # Load processed AnnData object
        print("Loading AnnData object...")
        adata_pp = sc.read_h5ad(os.path.join(working_dir, 'adata_pp.h5ad'))
        
        # Verify required fields
        required_columns = ['final_annotation', 'run_accession', 'series_id']
        for col in required_columns:
            if col not in adata_pp.obs.columns:
                raise ValueError(f"Missing column: {col}")
        
        # Process barcodes
        adata_pp.obs['original_barcode'] = adata_pp.obs_names.str.split('-').str[0]
        
        
        #==========================================#
        #  PHASE 1: Prepare inputs and run SplitBam #
        #==========================================#
        if RESTART_FROM_STEP <= 1:
            print("\n" + "="*80)
            print("PHASE 1: Preparing inputs and running SplitBam")
            print("="*80)
            
            run_accessions = adata_pp.obs['run_accession'].unique()
            scomatic_args = []
            
            for run_acc in tqdm(run_accessions, desc='Preparing samples'):
                # Setup directories
                mask = adata_pp.obs['run_accession'] == run_acc
                subset = adata_pp[mask, :]
                series_id = subset.obs['series_id'].iloc[0]
                
                cellranger_dir = os.path.join(working_dir, 'fastq', series_id, run_acc, f"{run_acc}_S1_L001_", "outs")
                original_bam = os.path.join(cellranger_dir, "possorted_genome_bam.bam")
                
                # Skip if files missing
                if not all(os.path.exists(p) for p in [cellranger_dir, original_bam]):
                    print(f"Skipping {run_acc} - files missing")
                    continue
                
                # Create filtered BAM
                scomatic_dir = os.path.join(cellranger_dir, 'SComatic')
                os.makedirs(scomatic_dir, exist_ok=True)
                
                filtered_bam = os.path.join(scomatic_dir, f"filtered_{run_acc}.bam")
                valid_barcodes = set(subset.obs['original_barcode'])
                
                # Skip if already processed
                if not os.path.exists(filtered_bam):
                    if not create_filtered_bam(original_bam, filtered_bam, valid_barcodes):
                        continue
                
                # Create annotation file
                annotation_file = os.path.join(scomatic_dir, 'cell_barcode_annotation.tsv')
                if not os.path.exists(annotation_file):
                    annotation_df = pd.DataFrame({
                        'Index': subset.obs['original_barcode'],
                        'Cell_type': subset.obs['final_annotation'].str.replace(r'[~.`!@#$%^&*(){|}/\\:;"\'<>?,=+\s]', '_', regex=True)
                    })
                    annotation_df.to_csv(annotation_file, sep='\t', index=False)
                
                scomatic_args.append((filtered_bam, annotation_file, run_acc, scomatic_dir, scomatic_scripts_dir))
    
            # Run SplitBam in parallel
            print(f"\nRunning SplitBam for {len(scomatic_args)} samples...")
            with ThreadPool(min(NCPUS, len(scomatic_args))) as pool:
                split_results = list(tqdm(
                    pool.imap(run_scomatic, scomatic_args),
                    total=len(scomatic_args),
                    desc="Running SplitBam"
                ))
            
            # Identify valid samples
            valid_samples = [args[2] for args, success in zip(scomatic_args, split_results) if success]
            print(f"SplitBam completed: {len(valid_samples)}/{len(scomatic_args)} successful")
        
        #==========================================#
        #  PHASE 2: BaseCellCounter Processing      #
        #==========================================#
        if RESTART_FROM_STEP <= 2:
            print("\n" + "="*80)
            print("PHASE 2: BaseCellCounter Processing")
            print("="*80)
            
            basecell_counter_args = []
            for run_acc in valid_samples:
                mask = adata_pp.obs['run_accession'] == run_acc
                subset = adata_pp[mask, :]
                series_id = subset.obs['series_id'].iloc[0]
                
                scomatic_dir = os.path.join(working_dir, 'fastq', series_id, run_acc, f"{run_acc}_S1_L001_", "outs", "SComatic")
                
                # Find split BAM files
                bam_files = glob.glob(os.path.join(scomatic_dir, f"{run_acc}.*.bam"))
                if not bam_files:
                    bam_files = [f for f in glob.glob(os.path.join(scomatic_dir, "*.bam")) 
                                if not f.endswith(f"filtered_{run_acc}.bam")]
                if not bam_files:
                    print(f"No split BAMs found for {run_acc}")
                    continue
                    
                # Create output directory
                basecell_counts_dir = os.path.join(scomatic_dir, 'BaseCellCounts')
                os.makedirs(basecell_counts_dir, exist_ok=True)
                
                for bam in bam_files:
                    basecell_counter_args.append((bam, ref_genome, basecell_counts_dir, scomatic_scripts_dir))
            
            # Run BaseCellCounter
            if basecell_counter_args:
                basecell_results = run_basecell_counter_parallel(basecell_counter_args, NCPUS)
            else:
                print("No BAM files to process")
                basecell_results = []
            
            # Save checkpoint
            checkpoint_path = os.path.join(working_dir, "basecell_counter_args.pkl")
            with open(checkpoint_path, 'wb') as f:
                pickle.dump((basecell_counter_args, basecell_results), f)
        
        #==========================================#
        #  PHASE 3: Merge Counts                   #
        #==========================================#
        if RESTART_FROM_STEP <= 3:
            print("\n" + "="*80)
            print("PHASE 3: Merging Counts")
            print("="*80)
            
            merge_results = []
            for run_acc in valid_samples:
                mask = adata_pp.obs['run_accession'] == run_acc
                subset = adata_pp[mask, :]
                series_id = subset.obs['series_id'].iloc[0]
                
                scomatic_dir = os.path.join(working_dir, 'fastq', series_id, run_acc, f"{run_acc}_S1_L001_", "outs", "SComatic")
                basecell_counts_dir = os.path.join(scomatic_dir, 'BaseCellCounts')
                
                if not os.path.exists(basecell_counts_dir):
                    print(f"Skipping {run_acc} - BaseCellCounts dir missing")
                    continue
                    
                merged_counts_dir = os.path.join(scomatic_dir, 'MergedCounts')
                os.makedirs(merged_counts_dir, exist_ok=True)
                
                output_file = os.path.join(merged_counts_dir, f"{run_acc}.BaseCellCounts.AllCellTypes.tsv")
                success = run_merge_counts(basecell_counts_dir, output_file, scomatic_scripts_dir, run_acc)
                merge_results.append(success)
            
            print(f"Merged counts completed: {sum(merge_results)}/{len(merge_results)} successful")
        
        #==========================================#
        #  PHASE 4: Variant Calling                #
        #==========================================#
        if RESTART_FROM_STEP <= 4:
            print("\n" + "="*80)
            print("PHASE 4: Variant Calling")
            print("="*80)
            
            variant_step1_args = []
            variant_step2_args = []
            variant_valid_samples = []
            
            for run_acc in valid_samples:
                mask = adata_pp.obs['run_accession'] == run_acc
                subset = adata_pp[mask, :]
                series_id = subset.obs['series_id'].iloc[0]
                
                scomatic_dir = os.path.join(working_dir, 'fastq', series_id, run_acc, f"{run_acc}_S1_L001_", "outs", "SComatic")
                merged_counts_file = os.path.join(scomatic_dir, 'MergedCounts', f"{run_acc}.BaseCellCounts.AllCellTypes.tsv")
                
                if not os.path.exists(merged_counts_file) or os.path.getsize(merged_counts_file) == 0:
                    print(f"Skipping {run_acc} - merged counts missing")
                    continue
                    
                variant_calling_dir = os.path.join(scomatic_dir, 'VariantCalling')
                os.makedirs(variant_calling_dir, exist_ok=True)
                
                output_prefix_step1 = os.path.join(variant_calling_dir, run_acc)
                variant_step1_args.append((merged_counts_file, output_prefix_step1, ref_genome, scomatic_scripts_dir))
                variant_step2_args.append((f"{output_prefix_step1}.calling.step1.tsv", output_prefix_step1, editing_sites, pon_file, scomatic_scripts_dir))
                variant_valid_samples.append(run_acc)
            
            # Process variant calling in parallel
            variant_results = []
            with ThreadPool(min(NCPUS, len(variant_step1_args))) as pool:
                # Step 1
                step1_results = list(tqdm(
                    pool.imap(run_variant_calling_step1, variant_step1_args),
                    total=len(variant_step1_args),
                    desc="Variant Calling Step 1"
                ))
                
                # Step 2 (only for successful step1)
                step2_args = [args for args, success in zip(variant_step2_args, step1_results) if success]
                step2_results = list(tqdm(
                    pool.imap(run_variant_calling_step2, step2_args),
                    total=len(step2_args),
                    desc="Variant Calling Step 2"
                ))
                
                # Combine results
                variant_results = [s1 and s2 for s1, s2 in zip(step1_results, step2_results)]
            
            # Update valid samples
            final_valid_samples = [run_acc for run_acc, success in zip(variant_valid_samples, variant_results) if success]
            print(f"Variant calling completed: {len(final_valid_samples)}/{len(variant_valid_samples)} successful")
        
        #==========================================#
        #  PHASE 5: BED Filtering                  #
        #==========================================#
        if RESTART_FROM_STEP <= 5:
            print("\n" + "="*80)
            print("PHASE 5: BED Filtering")
            print("="*80)
            
            bed_filter_args = []
            for run_acc in final_valid_samples:
                mask = adata_pp.obs['run_accession'] == run_acc
                subset = adata_pp[mask, :]
                series_id = subset.obs['series_id'].iloc[0]
                
                scomatic_dir = os.path.join(working_dir, 'fastq', series_id, run_acc, f"{run_acc}_S1_L001_", "outs", "SComatic")
                variant_calling_dir = os.path.join(scomatic_dir, 'VariantCalling')
                filtered_dir = os.path.join(scomatic_dir, 'FilteredVariants')
                os.makedirs(filtered_dir, exist_ok=True)
                
                # Find step2 files
                step2_files = glob.glob(os.path.join(variant_calling_dir, "*.calling.step2.tsv"))
                for variant_file in step2_files:
                    filename = os.path.basename(variant_file)
                    filtered_file = os.path.join(filtered_dir, filename.replace('.step2.tsv', '.filtered.tsv'))
                    bed_filter_args.append((variant_file, bed_file, filtered_file))
            
            # Run BED filtering
            bed_results = []
            with ThreadPool(min(NCPUS, len(bed_filter_args))) as pool:
                bed_results = list(tqdm(
                    pool.starmap(filter_with_bed, bed_filter_args),
                    total=len(bed_filter_args),
                    desc="BED Filtering"
                ))
            
            print(f"BED filtering completed: {sum(bed_results)}/{len(bed_results)} successful")
        
        #==========================================#
        #  PHASE 6: Final Processing               #
        #==========================================#
        if RESTART_FROM_STEP <= 6:
            # Ensure valid_samples exists for consistency
            if 'valid_samples' not in locals():
                valid_samples = final_valid_samples if final_valid_samples else []
            
            # Pre-flight validation checks
            print("\n" + "="*80)
            print("PRE-FLIGHT CHECKS FOR PHASE 6")
            print("="*80)
            
            # Check 1: Validate SingleCellGenotype.py has been modified
            singlecell_script = os.path.join(scomatic_scripts_dir, 'SingleCellGenotype', 'SingleCellGenotype.py')
            if os.path.exists(singlecell_script):
                with open(singlecell_script, 'r') as f:
                    script_content = f.read()
                    if 'Total_depth' not in script_content:
                        print("WARNING: SingleCellGenotype.py may not contain Total_depth modifications!")
                        print(f"  Location: {singlecell_script}")
                        print("  Expected: Script should contain 'Total_depth' in header and output")
                        response = input("\nContinue anyway? (yes/no): ").lower()
                        if response != 'yes':
                            print("Aborted. Please update SingleCellGenotype.py and try again.")
                            return
                    else:
                        print("Check 1: SingleCellGenotype.py contains Total_depth modifications")
            else:
                print(f"WARNING: Could not verify SingleCellGenotype.py at {singlecell_script}")
            
            # Check 2: Verify we have samples to process
            if not final_valid_samples or len(final_valid_samples) == 0:
                print("ERROR: No valid samples available for Phase 6")
                print("Previous phases must complete successfully first.")
                return
            
            print(f"Check 2: Processing {len(final_valid_samples)} samples")
            print("="*80)
            
            print("\n" + "="*80)
            print("PHASE 6: Final Processing")
            print("="*80)
            
            # Process samples in parallel
            if not all_sample_args:
                print("Preparing sample arguments for Phase 6...")
                all_sample_args = [
                    prepare_sample_args(
                        run_acc=run_acc,
                        adata_pp=adata_pp,
                        working_dir=working_dir,
                        scomatic_scripts_dir=scomatic_scripts_dir,
                        ref_genome=ref_genome
                    ) for run_acc in final_valid_samples
                ]
                all_sample_args = [args for args in all_sample_args if args is not None]
                print(f"After validation, {len(all_sample_args)} samples ready for Phase 6 processing")
                
                if len(all_sample_args) == 0:
                    print("ERROR: No samples passed validation for Phase 6")
                    print("Check that Phase 5 completed successfully for all samples")
                    return
            
            # Run Phase 6 processing
            final_results = []
            with ThreadPool(min(NCPUS, len(all_sample_args))) as pool:
                final_results = list(tqdm(
                    pool.imap(process_sample, all_sample_args),
                    total=len(all_sample_args),
                    desc="Phase 6 Processing"
                ))
            
            # Phase 6 diagnostics
            print("\n" + "="*80)
            print("PHASE 6 COMPLETION DIAGNOSTICS")
            print("="*80)
            
            successful_samples = 0
            failed_samples = 0
            
            for result in final_results:
                run_acc = result.get('run_acc', 'unknown')
                success = result.get('success', False)
                
                if success:
                    successful_samples += 1
                else:
                    failed_samples += 1
                    print(f"\nSample {run_acc}: FAILED")
                    print(f"  Callable sites: {result.get('callable_sites', 'N/A')}")
                    print(f"  Sites per cell: {result.get('sites_per_cell', 'N/A')}")
                    print(f"  Genotypes: {result.get('single_cell_genotypes', 'N/A')}")
                    print(f"  Filtered genotypes: {result.get('filtered_genotypes', 'N/A')}")
                    
                    if 'error' in result:
                        print(f"  Error: {result['error']}")
            
            print("\n  Phase 6 Summary:")
            print(f"  Successful: {successful_samples}/{len(final_results)}")
            print(f"  Failed: {failed_samples}/{len(final_results)}")
            
            # Validate genotype file structure
            print("\n" + "="*80)
            print("VALIDATING GENOTYPE FILE STRUCTURE")
            print("="*80)
            
            validation_passed = True
            for sample_args in all_sample_args:
                run_acc = sample_args['run_acc']
                single_cell_dir = sample_args['directories']['single_cell']
                
                if os.path.exists(single_cell_dir):
                    genotype_files = glob.glob(os.path.join(single_cell_dir, "*.single_cell_genotype.tsv"))
                    
                    if genotype_files:
                        # Check first file only
                        test_file = genotype_files[0]
                        try:
                            test_df = pd.read_csv(test_file, sep='\t', nrows=5)
                            expected_cols = ['#CHROM', 'Start', 'End', 'REF', 'ALT_expected', 
                                           'Cell_type_expected', 'Num_cells_expected', 'CB', 
                                           'Cell_type_observed', 'Base_observed', 'Num_reads', 'Total_depth']
                            
                            missing_cols = [col for col in expected_cols if col not in test_df.columns]
                            
                            if missing_cols:
                                print(f"ERROR: {run_acc} - Missing columns: {missing_cols}")
                                print(f"  File: {test_file}")
                                print(f"  Found columns: {test_df.columns.tolist()}")
                                validation_passed = False
                            else:
                                print(f"Validated: {run_acc} - All required columns present including Total_depth")
                                
                        except Exception as e:
                            print(f"ERROR reading {test_file}: {e}")
                            validation_passed = False
                    else:
                        print(f"WARNING: {run_acc} - No genotype files found")
            
            if not validation_passed:
                print("\nWARNING: Some genotype files are missing required columns")
                print("The filtering step may fail for these samples")
            
            # Compute trinucleotide context background
            print("\n" + "="*80)
            print("Computing trinucleotide context background...")
            print("="*80)
            
            trinuc_step1_files = []
            for sample_args in all_sample_args:
                step1_file = sample_args['step1_file']
                if os.path.exists(step1_file) and os.path.getsize(step1_file) > 0:
                    trinuc_step1_files.append(step1_file)
            
            if trinuc_step1_files:
                output_trinuc = os.path.join(working_dir, 'trinucleotide_background.tsv')
                trinuc_args = (trinuc_step1_files, output_trinuc, scomatic_scripts_dir)
                success = run_trinucleotide_context(trinuc_args)
                if success:
                    print(f"Trinucleotide background saved to {output_trinuc}")
                else:
                    print("Failed to compute trinucleotide background")
            else:
                print("No valid step1 files found for trinucleotide context")
            
            # Combine mutation data
            print("\n" + "="*80)
            print("Combining mutation data from all samples...")
            print("="*80)
            
            all_samples_data = []
            for result in final_results:
                if result.get('combined_data') is not None and not result['combined_data'].empty:
                    all_samples_data.append(result['combined_data'])
            
            print(f"Collected data from {len(all_samples_data)} samples")
            
            if all_samples_data:
                # Show preview of combined data
                for i, df in enumerate(all_samples_data[:3]):  # Show first 3 samples
                    print(f"\nSample {i+1} preview:")
                    print(f"  Rows: {df.shape[0]:,} variants")
                    print(f"  Unique cells: {df['CB'].nunique():,}")
                    print(f"  Columns: {df.columns.tolist()}")
                    if len(df) > 0:
                        print(f"  Trinuc example: REF_TRI={df['REF_TRI'].iloc[0]}, ALT_TRI={df['ALT_TRI'].iloc[0]}")
                        if 'Total_depth' in df.columns:
                            print(f"  Total_depth range: {df['Total_depth'].min()}-{df['Total_depth'].max()}")
                        else:
                            print("  WARNING: Total_depth column missing!")
                
                if len(all_samples_data) > 3:
                    print(f"\n... and {len(all_samples_data) - 3} more samples")
                
                # Combine all samples
                final_combined = pd.concat(all_samples_data, ignore_index=True)
                final_output_path = os.path.join(working_dir, 'all_samples.single_cell_genotype.filtered.tsv')
                final_combined.to_csv(final_output_path, sep='\t', index=False)
                
                print("\nCombined dataset statistics:")
                print(f"  Total variants: {len(final_combined):,}")
                print(f"  Unique cells: {final_combined['CB'].nunique():,}")
                print(f"  Saved to: {final_output_path}")
            else:
                print("\nWARNING: No mutation data collected from any samples")
                print("Check Phase 6 diagnostics above for failures")
        
        #==========================================#
        #  PHASE 7: Callable Sites                 #
        #==========================================#
        if RESTART_FROM_STEP <= 7:
            if final_valid_samples is None or len(final_valid_samples) == 0:
                print("\nERROR: No valid samples available for Phase 7")
                print("Phase 6 must complete successfully before running Phase 7")
            else:
                print("\n" + "="*80)
                print("PHASE 7: Generating Callable Sites")
                print("="*80)
                
                cell_annotations_path = os.path.join(working_dir, 'cell_annotations.txt')
                
                if not os.path.exists(cell_annotations_path):
                    print(f"ERROR: Cell annotations file not found: {cell_annotations_path}")
                else:
                    cell_annotations = pd.read_csv(cell_annotations_path, sep='\t')
                    
                    callable_path = generate_complete_callable_sites(
                        working_dir, final_valid_samples, adata_pp, cell_annotations
                    )
                    
                    if callable_path and os.path.exists(callable_path):
                        print(f"Callable sites generated successfully: {callable_path}")
                    else:
                        print("WARNING: Callable sites generation may have failed")
        
        #==========================================#
        #  PHASE 8: Final Verification             #
        #==========================================#
        if RESTART_FROM_STEP <= 8:
            print("\n" + "="*80)
            print("PHASE 8: Final Verification")
            print("="*80)
            
            # Define expected output paths
            final_output_path = os.path.join(working_dir, 'all_samples.single_cell_genotype.filtered.tsv')
            callable_path = os.path.join(working_dir, 'CombinedCallableSites', 'complete_callable_sites.tsv')
            output_trinuc = os.path.join(working_dir, 'trinucleotide_background.tsv')
            cell_annotations_path = os.path.join(working_dir, 'cell_annotations.txt')
            
            # Check which files exist
            files_status = {
                'mutations': os.path.exists(final_output_path),
                'callable_sites': os.path.exists(callable_path),
                'trinucleotide': os.path.exists(output_trinuc),
                'cell_annotations': os.path.exists(cell_annotations_path)
            }
            
            print("Output file status:")
            for file_type, exists in files_status.items():
                status = "EXISTS" if exists else "MISSING"
                print(f"  {file_type}: {status}")
            
            # Run verification if we have the required files
            if files_status['mutations'] and files_status['cell_annotations']:
                verify_pipeline_output(
                    final_output_path,
                    callable_path if files_status['callable_sites'] else "",
                    cell_annotations_path,
                    output_trinuc if files_status['trinucleotide'] else None
                )
            else:
                print("\nWARNING: Cannot run verification - required files missing")
                if not files_status['mutations']:
                    print("  Missing: all_samples.single_cell_genotype.filtered.tsv")
                if not files_status['cell_annotations']:
                    print("  Missing: cell_annotations.txt")
            
            # Final summary
            print("\n" + "="*80)
            print("PIPELINE COMPLETION SUMMARY")
            print("="*80)
            
            if all(files_status.values()):
                print("SUCCESS: All expected output files generated")
                print("\nKey outputs:")
                print(f"  1. Filtered mutations: {final_output_path}")
                print(f"  2. Callable sites: {callable_path}")
                print(f"  3. Trinucleotide background: {output_trinuc}")
                print("\nPipeline completed successfully!")
            else:
                print("PARTIAL SUCCESS: Some output files missing")
                missing = [k for k, v in files_status.items() if not v]
                print(f"  Missing files: {', '.join(missing)}")
                print("\nReview diagnostics above to identify issues")
        
    except Exception as e:
        print(f"\nCRITICAL ERROR in main pipeline: {str(e)}")
        import traceback
        traceback.print_exc()
    finally:
        # Stop resource monitoring
        monitor_resources.stop = True
        monitor_thread.join(timeout=10)


if __name__ == "__main__":
    main()
    
