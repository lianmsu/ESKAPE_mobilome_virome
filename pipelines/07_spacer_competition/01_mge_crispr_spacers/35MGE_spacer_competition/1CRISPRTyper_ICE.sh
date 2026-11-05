#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J CCtyper_ICE
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.01CCtyper_ICE-lpw.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 01_mge_crispr_spacers
# Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets)
# Manuscript      : Fig. 2B
# Original path   : lfat1t01:35MGE_spacer_competition/1CRISPRTyper_ICE.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate cctyper
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition

cctyper ./all_combined_ICEs.fa ./CCTyper_res_ICE -t 28 --prodigal single --no_grid
