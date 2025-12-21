#!/bin/bash

#SBATCH -J ADATA_SPATIAL                                              # Job name
#SBATCH -o /work/sdz852/WORKING/LOGS/ADATA_SPATIAL.o.log              # Name of the stdout output file
#SBATCH -e /work/sdz852/WORKING/LOGS/ADATA_SPATIAL.e.log              # Name of the stderr error file
#SBATCH --mail-user=jake.lehle@utsa.edu                                  # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 10-00:00:00                                                   # Time of job
#SBATCH -p compute2                                              # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 80                                                    # Total number of cores 80 max

# Startup scripts

module load anaconda3

conda activate slide-TCR-seq

conda run 
conda run -n slide-TCR-seq python /work/sdz852/WORKING/slide-TCR-seq/scripts/ADATA_PREP.py

exit
