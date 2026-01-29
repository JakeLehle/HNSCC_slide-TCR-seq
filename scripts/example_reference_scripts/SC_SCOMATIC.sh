#!/bin/bash

#SBATCH -J BRC_SCOMATIC                                              # Job name
#SBATCH -o /work/sdz852/WORKING/LOGS/BRC_MUT.o.log              # Name of the stdout output file
#SBATCH -e /work/sdz852/WORKING/LOGS/BRC_MUT.e.log              # Name of the stderr error file
#SBATCH --mail-user=jake.lehle@utsa.edu                                  # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 10-00:00:00                                                   # Time of job
#SBATCH -p compute2                                              # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 80                                                    # Total number of cores 80 max

# Load one of these
module load anaconda3

conda activate SComatic

cd /work/sdz852/WORKING/SC/fastq/Breast_cancer

conda run -n SComatic python /work/sdz852/WORKING/SC/fastq/Breast_cancer/SComatic_script.py

exit
