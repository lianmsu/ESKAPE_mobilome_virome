#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J bgc
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=job.01antiSMASH_lpw.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 03_resistome_and_bgc / 03_antismash
# Biosynthetic gene cluster prediction with antiSMASH
# Manuscript      : Fig. 1D
# Original path   : lfat1t01:33antismash/02antiSMASH.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/41antiSMASH
mkdir -p ./res1
source activate antiSMASH

antismash ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/41antiSMASH/merged_5455MAG.fa \
	--taxon bacteria \
	--cpus 40 \
	--databases ${HOME_ROOT:-/path/to/home}/software/miniforge3/envs/antiSMASH/lib/python3.11/site-packages/antismash/databases \
	--output-dir ./res1 \
	--allow-long-headers \
	--genefinding-tool prodigal-m

