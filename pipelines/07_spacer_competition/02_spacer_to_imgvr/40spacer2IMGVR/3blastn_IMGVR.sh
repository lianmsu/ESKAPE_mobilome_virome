#!/bin/bash
#SBATCH --job-name=blastnVR
#SBATCH --partition=fat
#SBATCH --nodes=1 
#SBATCH --get-user-env
#SBATCH --ntasks-per-node=48
#SBATCH --output=job-virus.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 02_spacer_to_imgvr
# BLAST of the host spacers against the IMG/VR database
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:40spacer2IMGVR/3blastn_IMGVR.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate blast2.16
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/40spacer2IMGVR

mkdir -p ./spacer2IMGVR

blastn -db ${HOME_ROOT:-/path/to/home}/06mobilome_predators/70_spacer2IMGVR/blastdb_IMGVR/IMGVR.db -query ./level4_minced_spacers.fa -out ./spacer2IMGVR/spacers_IMGVR.blastn.out -outfmt '6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send' -num_threads 48 -mt_mode 1 -evalue 1e-5 -word_size 8 -task blastn-short

cat ./spacer2IMGVR/spacers_IMGVR.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > spacers_IMGVR.blastn.1mis95cov.out
#cat spacers_virus.blastn.1mis95cov.out | cut -f2 | sort | uniq > targeted.virus.id


