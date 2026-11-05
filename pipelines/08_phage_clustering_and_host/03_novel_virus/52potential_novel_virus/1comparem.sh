
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 08_phage_clustering_and_host / 03_novel_virus
# Identification and comparison of the potential novel viruses
# Manuscript      : Fig. 5-6
# Original path   : lfat1t01:52potential_novel_virus/1comparem.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus

seqkit split -i ./39_phage.fasta -O ./split_39_phages

source activate comparem
comparem aai_wf --evalue 1e-5 --per_identity 30 --cpus 10 ./split_39_phages 39_phage_aai_output --file_ext fasta

