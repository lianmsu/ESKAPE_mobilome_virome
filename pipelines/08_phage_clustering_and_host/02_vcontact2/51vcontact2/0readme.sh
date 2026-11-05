
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 08_phage_clustering_and_host / 02_vcontact2
# vContact2 clustering for the companion comparison
# Manuscript      : Fig. 5
# Original path   : lfat1t01:51vcontact2/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

#NCBI refseq211 + public + ours，构建病毒网络
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/51vcontact2


cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus
source activate prokka

res_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus

contigs=${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus/ESKAPE_complete_phage.fa
echo "Running Prodigal"
prodigal -p meta -i "${contigs}" \
		-f gff -o "${res_path}/ESKAPE_complete_phage.gff" \
		-m -d "${res_path}/ESKAPE_complete_phage.fna" \
		-a "${res_path}/ESKAPE_complete_phage.faa"


cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/51vcontact2
source activate vContact2

cp ${res_path}/ESKAPE_complete_phage.faa ./

vcontact2_gene2genome \
--proteins ./ESKAPE_complete_phage.faa \
--output self_prepared_phage.csv \
--source-type Prodigal-FAA

vcontact2 --raw-proteins ./ESKAPE_complete_phage.faa --proteins-fp ./self_prepared_phage.csv --db 'ProkaryoticViralRefSeq211-Merged' --output-dir vConTACT2_Results
