# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 05_virus_discovery / 04_final_catalogue
# Filtering rules that define the final non-redundant ESKAPE virome
# Manuscript      : Fig. 3
# Original path   : lfat1t01:46final_virus/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/45virus_genomad/res/nr_phage_rep_seq_summary/nr_phage_rep_seq_virus.fna ./ESKAPE_all_phage.fa

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/45virus_genomad/res/nr_phage_rep_seq_summary/nr_phage_rep_seq_virus_summary.tsv ./ESKAPE_all_phage_genomad_summary.tsv

# 修正的一行命令版本
grep -i "DTR\|ITR" ./ESKAPE_all_phage_genomad_summary.tsv | cut -f1 | \
  seqkit grep -f - ./ESKAPE_all_phage.fa -o ESKAPE_complete_phage.fa
