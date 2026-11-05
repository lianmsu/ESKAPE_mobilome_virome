# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 05_virus_discovery / 02_checkv
# CheckV quality control of the geNomad virus candidates
# Manuscript      : Fig. 3A
# Original path   : lfat1t01:42genomad_checkv/02readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/42genomad_checkv

cat proviruses.fna viruses.fna > ESKAPE_genome_prophage.fa

# 排除序列名中含有 "plasmid" 的序列（不区分大小写）
seqkit seq -n -i ESKAPE_genome_prophage.fa | grep -vi plasmid | seqkit grep -f - ESKAPE_genome_prophage.fa -o ESKAPE_genome_prophage_filtered.fa

