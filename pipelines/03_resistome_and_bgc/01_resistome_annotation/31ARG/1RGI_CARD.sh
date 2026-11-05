#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J rgicard
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --output=job.02RGI_CARD_lpw.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 03_resistome_and_bgc / 01_resistome_annotation
# ARG annotation: RGI/CARD and diamond against SARG
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:31ARG/1RGI_CARD.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG
source activate RGI_CARD
mkdir -p ./0ARG_RGI

#rgi load --card_json ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/42ARG/database/ARG/card-data/card.json --local

rgi main --input_sequence ./4538genomes_prodigal.faa \
	--output_file ./0ARG_RGI/rgiargs \
	--local \
	--clean \
	-t protein \
	--num_threads 50

