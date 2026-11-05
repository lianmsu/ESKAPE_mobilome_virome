#!/bin/bash
#SBATCH --job-name=blastnVIRE
#SBATCH --partition=fat
#SBATCH --nodes=1 
#SBATCH --get-user-env
#SBATCH --ntasks-per-node=48
#SBATCH --output=job-vire.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 03_spacer_to_vire
# BLAST of the host spacers against the VIRE database
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:41spacer2VIRE/3blastn_VIRE.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate blast2.16
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/41spacer2VIRE

mkdir -p ./spacer2VIRE

blastn -db ${PROJECT_ROOT:-/path/to/08ESKAPE}/41spacer2VIRE/blastdb_VIRE/VIRE.db -query ./level4_minced_spacers.fa -out ./spacer2VIRE/spacers_VIRE.blastn.out -outfmt '6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send' -num_threads 48 -mt_mode 1 -evalue 1e-5 -word_size 8 -task blastn-short

cat ./spacer2VIRE/spacers_VIRE.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > spacers_VIRE.blastn.1mis95cov.out
#cat spacers_virus.blastn.1mis95cov.out | cut -f2 | sort | uniq > targeted.virus.id


