# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 02_spacer_to_imgvr
# BLAST of the host spacers against the IMG/VR database
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:40spacer2IMGVR/2makeblastdb.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/40spacer2IMGVR

ln -s ${HOME_ROOT:-/path/to/home}/06mobilome_predators/70_spacer2IMGVR/blastdb_IMGPR/

# #运行以下命令查找重复的 ID：
# grep "^>" IMGVR_all_nucleotides.fna | sort | uniq -d
# #使用 seqkit 工具自动去重（保留第一条出现的序列）：
# seqkit rmdup -s IMGVR_all_nucleotides.fna -o IMGVR_all_nucleotides_dedup.fna

# #构建成功
# nohup makeblastdb -in ${HOME_ROOT:-/path/to/home}/06mobilome_predators/2IMG_VR_4/IMGVR_all_nucleotides.fna -dbtype nucl -out ./blastdb_IMGVR/IMGVR.db -max_file_sz '4GB' > nohup_makeblastdb_IMGVR.log 2>&1 &

# #报错序列id太长
# nohup makeblastdb -in ${HOME_ROOT:-/path/to/home}/06mobilome_predators/2IMG_VR_4/IMGVR_all_nucleotides.fna -dbtype nucl -parse_seqids -out ./blastdb_IMGVR/IMGVR.db -max_file_sz '4GB' > nohup_makeblastdb_IMGVR.log 2>&1 &

# #报错存在重复序列，
# nohup makeblastdb -in ${HOME_ROOT:-/path/to/home}/06mobilome_predators/2IMG_VR_4/IMGVR_all_nucleotides.fna -dbtype nucl -out ./blastdb_IMGVR/IMGVR.db -hash_index -max_file_sz '4GB' > nohup_makeblastdb_IMGVR.log 2>&1 &

# #sz服务器无法makeblastdb建库，在bb建好copy过来的。

# mkdir -p ./blastdb_IMGPR

# nohup makeblastdb -in ${HOME_ROOT:-/path/to/home}/06mobilome_predators/2IMG_PR/IMGPR_nucl.fna -dbtype nucl -out ./blastdb_IMGPR/IMGPR.db -hash_index -max_file_sz '4GB' > nohup_makeblastdb_IMGPR.log 2>&1 &

