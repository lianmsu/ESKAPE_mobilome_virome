# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 05_virus_discovery / 03_clustering
# Combine the virus candidates from all sources
# Manuscript      : Fig. 3
# Original path   : lfat1t01:43combined_all_virus/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/43combined_all_virus

#IMGVR
cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/40spacer2IMGVR/spacers_IMGVR.blastn.1mis100cov.fa ./
#VIRE
cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/41spacer2VIRE/spacers_VIRE.blastn.1mis100cov.fa ./
#vhdb
cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/vhdb_ESKAPE_virus.fasta ./vhdb_ESKAPE_virus.fa
#provirus
cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/42genomad_checkv/ESKAPE_genome_prophage_filtered.fa ./
#整理病毒序列名对应的宿主


#
cat ./*.fa > all_ESKAPE_phage.fa

