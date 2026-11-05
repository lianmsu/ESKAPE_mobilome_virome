#!/bin/bash
#SBATCH -o job.checkm2.%j.%N.out 
#SBATCH --partition=cpu
#SBATCH -J checkm2
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 01_genome_catalogue / 02_quality_and_taxonomy
# Genome quality with CheckM2
# Manuscript      : Fig. 1B
# Original path   : lfat1t01:1genome/3checkm2/checkm2.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate checkm2
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/3checkm2

checkm2 predict --threads 25 -i ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/renamed_4538genomes -o ./res -x .fa --database_path ${DB_ROOT:-/path/to/databases}/CheckM2_database/uniref100.KO.1.dmnd

#checkm2 predict --threads 10  -i ${HOME_ROOT:-/path/to/home}/data/2269_gw_ar_bins -o ${HOME_ROOT:-/path/to/home}/03GW_archaea/02checkm2 -x .fna
