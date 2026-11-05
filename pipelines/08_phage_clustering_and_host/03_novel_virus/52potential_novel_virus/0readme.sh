# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 08_phage_clustering_and_host / 03_novel_virus
# Identification and comparison of the potential novel viruses
# Manuscript      : Fig. 5-6
# Original path   : lfat1t01:52potential_novel_virus/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus

cat 1_potential_novel_virus.txt | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus/ESKAPE_complete_phage.fa > ./potential_novel_virus.fa

./ViPTreeGen ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/novel_virus_with_vhdb.fa ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/try

#使用NCBI Entrez得到密切相关的64个病毒
cat ./potential_novel_virus.fa ./64_phage.fasta > ./71_phage.fasta

cat ./potential_novel_virus.fa ./32_phage.fasta > ./39_phage.fasta
