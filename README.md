# HNSCC_slide-TCR-seq
This is a repository which organizes all scripts involved in processing HNSCC slide-TCR-seq data

# 2025-12-21
Initial upload of scripts I put together. This covers making a STAR reference genome and then trying to demultiplex the spatial data, starting with raw .blc files to fastqs and then pushing them through STARsolo to make up the spatial mapped count matrix. 

RCLONE.sh 
(moves the files over from dropbox)
MAKE_CONDA_ENV.sh 
(gets the conda env setup)
environment_slide-TCR-seq.yml 
(conda env)
REFERNCE_GENOME.sh 
(downloads GRCh38 and makes a STAR reference genome)
ALIGN_READS.sh
(Aligns the spatial reads using STARsolo to make up the spatially mapped count matrix. NOTE: The script works but its missing ~95% of the barcodes that should be in the whitelist. Looks like I don’t have the TCR seq bcl data for H52J2DMXY.)
 
I hit some snags looks like some of the files for the TCR sequences needed to be moved over in the repo to finish getting all this set up but I had access to the fully processed data. I can see by looking at some of the processed samples it looks like the Liu Lab used Drop-seq Tools so it’s possible I might need to switch over to that pipeline but I should probably talk with Mehdi Borji to get a better understanding of what all he has set up as I go along. I can see he has some repos where he does reads in spatial data but they looked more custom for his stuff and might take a bit to unpack.

Working with the processed data I went through and did the standard QC for working with single cell and removed low quality cels (low number of reads or unique genes, about 10K of the ~45K cells, leaving about ~35K cells per sample to work with). I’ll review and update the QC figures to dial them in cause the formatting was off. 

ADATA_PREP.sh/py (Makes the adata objects for each puck and then does general QC from the theislabsingle-cell-best-practices that I use from my other pipelines)

