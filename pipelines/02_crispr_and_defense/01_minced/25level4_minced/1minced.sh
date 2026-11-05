#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J mc
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.01minced_lpw.%j.%N.out
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 01_minced
# minced spacer extraction on the level-4 CRISPR arrays
# Manuscript      : Fig. 2B
# Original path   : lfat1t01:25level4_minced/1minced.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/25level4_minced
mkdir -p ./minced_res
source activate minced

minced -spacers ${PROJECT_ROOT:-/path/to/08ESKAPE}/24CRISPRCasTyper/all_CRISPR_level4.fa ./minced_res/level4_minced.crisprs ./minced_res/level4_minced.gff

