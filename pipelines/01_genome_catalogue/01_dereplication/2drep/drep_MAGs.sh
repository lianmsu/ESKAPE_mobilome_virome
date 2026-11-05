#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J drep2
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=500GB
#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 01_genome_catalogue / 01_dereplication
# dRep dereplication of the 11 947 ESKAPE genomes
# Manuscript      : Fig. 1A
# Original path   : pkusz_CPU:2drep/drep_MAGs.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate drep 
#module load CheckM/1.1.2

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep
#合并到一起
mkdir -p ./tmp

#ln -s ~/07deepwell/05bins/all_single_sample_bins/*.fa tmp
#ln -s ~/07deepwell/05bins/all_multiple_sample_bins/*.fa tmp

#去冗余
dRep dereplicate res2_ESKAPE -g ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_genomes/*.fa -sa 0.999 -nc 0.30 -p 45 -comp 50 -con 10
