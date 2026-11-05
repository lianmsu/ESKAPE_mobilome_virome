#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J finder
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --output=job.01CASfinder_lpw.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 02_crisprcasfinder
# CRISPR-Cas detection with CRISPRCasFinder (level 4 arrays)
# Manuscript      : Fig. 2B
# Original path   : lfat1t01:23CRISPRCasFinder/1CRISPRCasFinder.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/23CRISPRCasFinder
source activate crisprcasfinder
cd ./CRISPRCasFinder-master

perl CRISPRCasFinder.pl -in ../minced_crispr_contig.fa -levelMin 4 -q

#perl CRISPRCasFinder.pl -in ${HOME_ROOT:-/path/to/home}/06mobilome_predators/6_CRISPRCasFinder/minced_crispr_contig.fasta -levelMin 4 -cas -meta -q -cpuMacSyFinder 18
