# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 08_phage_clustering_and_host / 04_multi_host_phage
# Multi-host phages and their host ranges
# Manuscript      : Fig. 4
# Original path   : lfat1t01:53multi_host_phage/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/53multi_host_phage

#从本地R语言筛选出多宿主病毒的metadata，complete

cat ./En_Kp_EnKp_phage_complete.tsv | cut -f1 | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.fa > En_Kp_EnKp_phage_complete.fa

source activate prokka
res_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/53multi_host_phage
contigs=${PROJECT_ROOT:-/path/to/08ESKAPE}/53multi_host_phage/En_Kp_EnKp_phage_complete.fa
echo "Running Prodigal"
prodigal -p meta -i "${contigs}" \
		-f gff -o "${res_path}/En_Kp_EnKp_phage_complete.gff" \
		-m -d "${res_path}/En_Kp_EnKp_phage_complete.fna" \
		-a "${res_path}/En_Kp_EnKp_phage_complete.faa"

source activate vContact2
vcontact2_gene2genome \
--proteins ./En_Kp_EnKp_phage_complete.faa \
--output self_prepared_phage.csv \
--source-type Prodigal-FAA

vcontact2 --raw-proteins ./En_Kp_EnKp_phage_complete.faa --proteins-fp ./self_prepared_phage.csv --db 'ProkaryoticViralRefSeq211-Merged' --output-dir vConTACT2_Results

