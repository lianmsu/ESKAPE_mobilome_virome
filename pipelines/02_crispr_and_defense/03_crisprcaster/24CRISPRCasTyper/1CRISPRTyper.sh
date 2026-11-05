#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J CCtyper
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.01CCtyper_lpw.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 03_crisprcaster
# CRISPR-Cas subtyping with CRISPRCasTyper
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:24CRISPRCasTyper/1CRISPRTyper.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate cctyper
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/24CRISPRCasTyper

cctyper ./all_CRISPR_level4.fa ./CCTyper_res -t 28 --prodigal single --no_grid
