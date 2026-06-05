#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 10:58:12 2024

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
import pickle

with open(os.path.join(working_dir, 'metadata/dictionary_file_filtered_post_adata.pkl'), 'rb') as pkl_file:
     gse_dict = pickle.load(pkl_file)
     print('Dictionary loaded successfully')

#%% Define all the functions

import pysam
import regex as re
import collections.abc
from scipy.sparse import csr_matrix
from scipy.io import mmwrite
import csv
import logging
import os, sys
    
def extract_ids(bamfile, krakenfile): 
    """
    Builds a nested dictionary with KRAKEN2 taxonomy code for each transcript and the cell it belongs to.
    Input:  Output from KRAKEN2, .bam file with unmapped reads
    Output: {cellbarcode: {transcriptbarcode: krakentaxonomyID}}
    """
    line = 0
    skipped = 0
    # Store extracted information in nested dictionary {cellbarcode:{transcriptbarcode: taxonomyID}}
    nested_dict = {}
    
    # Iterate simultanously through bam and kraken file
    for sread,kread in zip(pysam.AlignmentFile(bamfile, "rb"),open(krakenfile,"r")):
        
        # count the total number of reads analysed
        line += 1
        
        # Check that read names in kraken and bam file match
        if sread.query_name != kread.split('\t')[1]:
            skipped += 1
            logging.warning("sam file read name and metagenomicsfile read name don't match and are therefore excluded: sam: {}, kraken: {}".format(sread.query_name, kread.split('\t')[1]))
            continue

        # Get cell barcode and UMI from bam file
        try:
            sread_CB = sread.get_tag('CB')
            sread_UB = sread.get_tag('UB')
        except:
            # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
            skipped += 1
            continue
            
        # Get taxonomy ID from kraken file
        kread_taxid = kread.split('\t')[2]
        if (type(kread_taxid) != int) and (kread_taxid.isdigit() == False):
            try:
                # sometimes, the taxonomy is name (taxid #), sometimes it's just the number
                kread_taxid = re.search('\(([^)]+)', kread_taxid).group(1)[6:]
            except:
                # in this case, something is wrong!
                logging.debug("Here is an error. TaxID: {}".format(kread_taxid))
                sys.exit()

        # Make nested dictionary with cells and transcripts
        if sread_CB in nested_dict:
            # If cell and transcript exist, add taxonomy ID to list
            if sread_UB in nested_dict[sread_CB]:
                nested_dict[sread_CB][sread_UB].append(kread_taxid)
            # Otherwise create transcript dictionary for cell
            else:
                nested_dict[sread_CB][sread_UB] = [kread_taxid]
        else:
            # if cell doesn't exist, create cell and transcript dictionary with kraken id
            nested_dict[sread_CB] = {sread_UB: [kread_taxid]}
    
    # Output control values
    logging.info("total reads: {}, skipped reads: {}".format(line,skipped))
    
    return nested_dict

def most_frequent(List):
    """Finds the most frequent element in a list"""
    return max(set(List), key = List.count)

def map_nested_dicts(ob, func):
    """ Applys a map to the inner item of nested dictionaries """
    for k, v in ob.items():
        if isinstance(v, collections.abc.Mapping):
            map_nested_dicts(v, func)
        else:
            ob[k] = func(v)

def twist_dict(nested):
    """ Make count dictionary with {cellbarcode : {taxonomyID : transcriptcount}} """
    newdict = {}
    for ckey, tdict in nested.items():
        for tkey, kvalue in tdict.items():
            if ckey in newdict:
                if kvalue in newdict[ckey]:
                    newdict[ckey][kvalue] += 1
                else:
                    newdict[ckey][kvalue] = 1
            else:
                newdict[ckey] = {kvalue: 1}
    return(newdict)

def dict2lists(nested):
    """ Returns lists for sparse matrix """
    rows = [] # cell coordinate
    columns = [] # taxonomy id coordinate
    values = [] # count

    cell_list = [] # same order as rows
    taxid_list = [] # same order as columns

    j = 0

    for ckey, taxdict in nested.items():
        for taxkey, count in taxdict.items():
            try:
                k = taxid_list.index(taxkey)
            except:
                taxid_list.append(taxkey)
                k = taxid_list.index(taxkey)
                
            rows.append(k)
            columns.append(j)
            values.append(count) 
            
        # increase cell coordinate by 1
        cell_list.append(ckey)
        j += 1
    
    return rows, columns, values, cell_list, taxid_list

def krakenID2dict(dbfile, taxid_list):
    """
    Get name for each taxonomy ID from kraken database
    """
    # iterate through inspect file and lookup taxonomy ids
    k=0
    taxdict = {'0': 'unclassified'}

    with open(dbfile) as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            # string to list
            line = line[:-1].split('\t')
            taxid_db = line[4]
            taxname = line[5].lstrip()
            
            if taxid_db in taxid_list:
                taxdict[taxid_db] = taxname
    
    return taxdict

def extract_taxref(file):
    """ 
    Extract taxonomy reference for each read.
    Input:  viral track output .bam file
    Output: dictionary with {readname: taxonomy ID}, list of unique taxonomy IDs
    """
    # extract taxref for each read
    tdict = {}
    line = 0
    skipped = 0
    taxref_list = set('0')
    
    for read in pysam.AlignmentFile(file, "rb"):
        # count the total number of reads analysed
        line += 1
        try:
            # Extract readname and taxonomy reference
            taxref = read.to_dict().get('ref_name').split('|')[1]
            taxref_list.add(taxref)
            tdict[read.query_name] = taxref
        except:
            # in case some reads are unmapped or don't work
            skipped += 1
    logging.info("Reads in ViralTrack output: {}, reads without taxonomy reference or that failed: {}".format(line, skipped))
    return(tdict, taxref_list)

def extract_bc(file):
    """ 
    Extracts cellbarcode and UMI for each readname
    Input:  unmapped .bam file
    Output: dictionary with {readname: [cellbarcode, UMI]}
    """
    # extract UB and CB for each read
    bcdict = {}
    line = 0
    skipped = 0

    for read in pysam.AlignmentFile(file, "rb"):
        # count the total number of reads analysed
        line += 1
        # Get cell barcode and UMI from bam file
        try:
            # Extract readname, cell barcode and UMI
            bcdict[read.query_name] = [read.get_tag('CB'),read.get_tag('UB')]
        except:
            # some reads don't have a cellbarcode or transcript barcode. They can be skipped.
            skipped += 1
            continue

    logging.info("Reads in original bam file: {}, reads without cellbarcode or UMI: {}".format(line, skipped))
    return(bcdict)


def mg2sc(bamfile, mgfile, dbfile, outdir):
    """ Main Function. 
    Creates a sparse matrix with transcript count per organism for each cell."""

    # Generate variables based on input
    try:
        os.mkdir(outdir + '/kraken2_filtered_feature_bc_matrix')
    except OSError as error:
        print(error)
        pass
    matrixfile = outdir + '/kraken2_filtered_feature_bc_matrix/matrix.mtx'
    cellfile = outdir + '/kraken2_filtered_feature_bc_matrix/barcodes.tsv'
    taxfile = outdir + '/kraken2_filtered_feature_bc_matrix/genes.tsv'
    dbfile = os.path.join(dbfile, 'inspect.txt')
    dbfile_out = outdir + '/kraken2_filtered_feature_bc_matrix/hierarchy.txt'

    # Extract taxonomy IDs for each transcript
    mg_dict = extract_ids(bamfile, mgfile)

    # Find most frequent taxonomy for each transcript
    map_nested_dicts(mg_dict, most_frequent)

    # Make sparse matrix
    rows, cols, vals, cell_list, taxid_list = dict2lists(twist_dict(mg_dict))
    sparsematrix = csr_matrix((vals, (rows, cols)))

    # Get ncbi name for taxonomy ID
    taxdict = krakenID2dict(dbfile, taxid_list)
    taxname_list = [taxdict[k] for k in taxid_list]

    # store sparse matrix
    mmwrite(matrixfile, sparsematrix)
    
    # Store list of cell barcodes
    with open(cellfile, 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\n')
        tsv_output.writerow(cell_list)
    
    # Store list of taxonomy IDs
    data = zip(taxid_list, taxname_list)
    with open(taxfile, 'w') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\t')
        for idx, tax  in data:
            tsv_output.writerow([idx, tax])
            
    # Store reference database hierarchy
    with open(dbfile) as f:
        with open(dbfile_out, "w") as f1:
            for line in f:
                f1.write(line) 


#%% Get the gse_dict

import subprocess
from os import cpu_count

NCPUS = cpu_count()

#For the human filtered dir
#k2_db=("kraken2/human_viral")

k2_db=("/work/sdz852/WORKING/SC/ALL_VIRUS/viral")

for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        
        bamfile_out = os.path.join(working_dir, "fastq", key, accession, accession+"_S1_L001_/outs", "possorted_genome_bam_unmapped.bam")
        out_dir = os.path.join(working_dir, "fastq", key, accession, accession+"_S1_L001_/outs")
        fqfile = os.path.join(working_dir, "fastq", key, accession, accession+"_S1_L001_/outs", "possorted_genome_bam_unmapped.fq")
        reportf = os.path.join(working_dir, "fastq", key, accession, accession+"_S1_L001_/outs", "possorted_genome_bam_krakenreport.txt")
        krakenoutfile = os.path.join(working_dir, "fastq", key, accession, accession+"_S1_L001_/outs", "possorted_genome_bam_output.kraken")
        
        cmd1 = "samtools view -@ " + str(NCPUS) + " -b -f 4 " + os.path.join(out_dir, "possorted_genome_bam.bam") \
            + " > " + bamfile_out
        os.makedirs(out_dir, exist_ok=True)
        proc1 = subprocess.Popen(cmd1, shell=True)
        proc1.wait()
        logging.info("Unmapped reads were extracted and saved to {}".format(bamfile_out))
        
        cmd2 = "samtools fastq -@ " + str(NCPUS) + " -n " \
            + bamfile_out \
            + " > " + fqfile
        proc2 = subprocess.Popen(cmd2, shell=True)
        proc2.wait()
        logging.info("FASTQ generated and saved to {}".format(fqfile))
        
        cmd3 = "kraken2 --threads " + str(NCPUS) + " --db " + os.path.join(os.environ['HOME'], k2_db) \
            + " --report " + reportf \
            + " " + fqfile \
            + " > " + krakenoutfile
        proc3 = subprocess.Popen(cmd3, shell=True)
        proc3.wait()
        logging.info("Kraken2 finished running, krakenoutput saved to {}".format(krakenoutfile))
        
        mg2sc(bamfile_out, krakenoutfile, \
              k2_db, \
              out_dir)
        logging.info("Sparse matrix with single cell information created")
        logging.info("Run finished.")
        


#%% gzip all the files so they look like the cell ranger output and save space
import gzip
import shutil
for key in gse_dict.keys():
    for accession in gse_dict[key]['run_accession']:
        out_dir = os.path.join(working_dir, "fastq", key, accession, accession+"_S1_L001_/outs/kraken2_filtered_feature_bc_matrix")
        files_list = ("matrix.mtx", "barcodes.tsv", "genes.tsv")
        for file in files_list:
            with open(os.path.join(out_dir, file), 'rb') as orginal_file:
                with gzip.open(os.path.join(out_dir, file+".gz"), 'wb') as zipped_file:
                    shutil.copyfileobj(orginal_file, zipped_file)
            os.remove(os.path.join(out_dir, file))
#%% End file
import sys

sys.exit()
