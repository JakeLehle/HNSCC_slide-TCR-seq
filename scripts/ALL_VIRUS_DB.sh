#!/bin/bash

#SBATCH -J BRC_K2                                              # Job name
#SBATCH -o /work/sdz852/WORKING/LOGS/BRC_K2.o.log              # Name of the stdout output file
#SBATCH -e /work/sdz852/WORKING/LOGS/BRC_K2.e.log              # Name of the stderr error file
#SBATCH --mail-user=jake.lehle@utsa.edu                                  # Send me an email when you are done
#SBATCH --mail-type=ALL
#SBATCH -t 10-00:00:00                                                   # Time of job
#SBATCH -p compute2                                              # Queue (partition) name
#SBATCH -N 1                                                     # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                                                     # Total # of mpi tasks (should be 1 for serial)
#SBATCH -c 80                                                    # Total number of cores 80 max


###### These sections have each been run but I'm leaving thme to make my life easy
#conda install bioconda::kraken2
#git clone https://github.com/DerrickWood/kraken2.git
#cd $HOME/kraken2
#wget https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
#tar -zxvf ncbi-blast-2.17.0+-x64-linux.tar.gz
#export PATH="/work/sdz852/WORKING/kraken2/ncbi-blast-2.16.0+/bin:$PATH"

# Startup scripts
module load anaconda3

#conda create -n kraken2_hu_vir
#conda activate kraken2_hu_vir
#conda install conda-forge::gxx_linux-64
#conda install bioconda::blast
#conda install conda-forge::glib


#Once the above has been run the first time you can skip it and start here
#conda activate kraken2_hu_vir

#cd ~/WORKING/kraken2

#./install_kraken2.sh /work/sdz852/WORKING/kraken2
#export PATH="/work/sdz852/WORKING/kraken2:$PATH"

###### Here is what you need to run to make the human viral database
#k2 download-taxonomy --db human_viral
#k2 add-to-library --threads 96 --db $HOME/WORKING/kraken2/human_viral --files '/master/jlehle/SC/ALL_VIRUS/ref/viral_ref_edited.fasta'
#k2 download-taxonomy --db /master/jlehle/WORKING/kraken2/human_viral
#k2 build --db human_viral --threads 96

export KRAKEN2_DB_PATH="/work/sdz852/WORKING/SC/ALL_VIRUS/"

# Then inspect the database by giving the name of the database for Kraken2 to look for
#kraken2-inspect --db human_viral > human_viral/inspect.txt

conda activate kraken2

# Finally run a python script to use the Kraken2 database on the single cell data you have collected
conda run -n kraken2 python3 /work/sdz852/WORKING/SC/fastq/Breast_cancer/Kraken2_script.py

exit

