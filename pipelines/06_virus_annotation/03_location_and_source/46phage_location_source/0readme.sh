# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 03_location_and_source
# Location and source metadata of the retained viruses
# Manuscript      : Fig. 4A
# Original path   : lfat1t01:46phage_location_source/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_location_source

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.tsv ./

ln -s ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE/vire_metadata_env_all_samples.tsv ./

ln -s ${HOME_ROOT:-/path/to/home}/06mobilome_predators/2IMG_VR_4/IMGVR_all_Sequence_information.tsv ./

head -n 1 ${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_location_source/ESKAPE_all_phage_info_filter.tsv > ./ESKAPE_IMGVR_phage.tsv
grep 'IMGVR' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_location_source/ESKAPE_all_phage_info_filter.tsv >> ./ESKAPE_IMGVR_phage.tsv
#在本地excel中进行处理

head -n 1 ${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_location_source/ESKAPE_all_phage_info_filter.tsv > ./ESKAPE_vire_phage.tsv
grep 'vire' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_location_source/ESKAPE_all_phage_info_filter.tsv >> ./ESKAPE_vire_phage.tsv
#提取文件vire_metadata_env_all_samples.tsv中第一列属于文件ESKAPE_vire_phage.tsv第一列的行
# 使用awk进行精确字段匹配
awk 'BEGIN {FS=OFS="\t"} 
     NR==FNR {ids[$1]; next} 
     FNR==1 || ($1 in ids)' \
     ESKAPE_vire_phage.tsv \
     vire_metadata_env_all_samples.tsv \
     > vire_matched_rows.tsv

tail -n+2 vire_matched_rows.tsv | cut -f2 > biosample_ids.txt

