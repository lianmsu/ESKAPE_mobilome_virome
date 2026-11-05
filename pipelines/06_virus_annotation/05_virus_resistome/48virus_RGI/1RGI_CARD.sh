#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J rgicard
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.02RGI_CARD_lpw.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 05_virus_resistome
# ARG annotation of the viruses with RGI/CARD
# Manuscript      : Fig. 3F
# Original path   : lfat1t01:48virus_RGI/1RGI_CARD.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/48virus_RGI

#cp -r ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG/localDB ./

source activate RGI_CARD
mkdir -p ./0ARG_RGI2

#rgi load --card_json ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/42ARG/database/ARG/card-data/card.json --local

rgi main --input_sequence ${PROJECT_ROOT:-/path/to/08ESKAPE}/47virus_defensefinder/res2/ESKAPE_all_phage_info_filter.prt \
	--output_file ./0ARG_RGI2/rgiargs \
	--local \
	--clean \
	-t protein \
	--num_threads 25

