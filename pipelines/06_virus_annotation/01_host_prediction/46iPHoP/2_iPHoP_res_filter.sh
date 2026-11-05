# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 01_host_prediction
# Host prediction with iPHoP and lifestyle with BACPHLIP
# Manuscript      : Fig. 3C-E
# Original path   : lfat1t01:46iPHoP/2_iPHoP_res_filter.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP

head -n1 ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv > ./iphop_phage2hostgenome.tsv

grep 'Enterococcus_B faecium' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | sort | uniq >> iphop_phage2hostgenome.tsv

grep 'Staphylococcus aureus' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | sort | uniq >> iphop_phage2hostgenome.tsv

grep 'Klebsiella pneumoniae' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | sort | uniq >> iphop_phage2hostgenome.tsv

grep 'Acinetobacter baumannii' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | sort | uniq >> iphop_phage2hostgenome.tsv

grep 'Pseudomonas aeruginosa' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | sort | uniq >> iphop_phage2hostgenome.tsv

grep 'g__Enterobacter' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb/Host_prediction_to_genome_m90.tsv | sort | uniq >> iphop_phage2hostgenome.tsv

#本地用R语言进行病毒序列的进一步筛选，得到ESKAPE_all_phage_info_filter.tsv
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP
cut -f1 ./ESKAPE_all_phage_info_filter.tsv | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus/checkv_clustered_genomad_checkv_filtered_phage.fa > ESKAPE_all_phage_info_filter.fa


