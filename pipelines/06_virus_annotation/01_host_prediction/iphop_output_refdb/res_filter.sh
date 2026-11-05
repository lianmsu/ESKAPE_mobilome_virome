# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 01_host_prediction
# Host prediction with iPHoP and lifestyle with BACPHLIP
# Manuscript      : Fig. 3C-E
# Original path   : lfat1t01:46iPHoP/iphop_output_refdb/res_filter.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb

grep 'Enterococcus_B faecium' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq | wc -l
#468

grep 'Staphylococcus aureus' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq | wc -l
#781

grep 'Klebsiella pneumoniae' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq | wc -l
#4773

grep 'Acinetobacter baumannii' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq | wc -l
#2727

grep 'Pseudomonas aeruginosa' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq | wc -l
#1622

grep 'g__Enterobacter' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq | wc -l
#4740

###############################################
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb
grep 'Enterococcus_B faecium' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq > phage_Ef.id

grep 'Staphylococcus aureus' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq > phage_Sa.id

grep 'Klebsiella pneumoniae' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq > phage_Kp.id

grep 'Acinetobacter baumannii' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq > phage_Ab.id

grep 'Pseudomonas aeruginosa' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq > phage_Pa.id

grep 'g__Enterobacter' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | cut -f1 | sort | uniq > phage_En.id

cat ./phage_*.id | sort | uniq | wc -l
cat ./phage_*.id | sort | uniq > ESKAPE_phage.id

