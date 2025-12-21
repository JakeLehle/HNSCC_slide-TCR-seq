#!/bin/bash

#SBATCH -J RCLONE_SLIDE-TCR-SEQ                                              # Job name
#SBATCH -o /work/sdz852/WORKING/LOGS/RCLONE_SPATIAL.o.log              # Name of the stdout output file
#SBATCH -e /work/sdz852/WORKING/LOGS/RCLONE_SPATIAL.e.log              # Name of the stderr error file
#SBATCH --mail-user=jake.lehle@utsa.edu                                  # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 10-00:00:00                                                   # Time of job
#SBATCH -p compute2                                              # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 80                                                    # Total number of cores 80 max

# Startup scripts

module load anaconda3

conda activate rclone

cd /work/sdz852/WORKING/slide-TCR-seq/

mkdir -p /work/sdz852/WORKING/slide-TCR-seq/fastq

rclone copy remote:/Jake\ Lehle/Spatial_Data_Sophia_Liu_HNSCC /work/sdz852/WORKING/slide-TCR-seq/fastq

exit
