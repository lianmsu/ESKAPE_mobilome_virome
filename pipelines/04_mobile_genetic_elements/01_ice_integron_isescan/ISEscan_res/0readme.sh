
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 01_ice_integron_isescan
# Insertion sequence detection with ISEScan
# Manuscript      : Fig. 2A
# Original path   : pkusz_CPU:6ISEscan/ISEscan_res/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

#合并路径${HOME_ROOT:-/path/to/home}/05plateau/23ISEscan_5455MAGs/results/5455MAG下所有.tsv结尾的文件到${HOME_ROOT:-/path/to/home}/05plateau/23ISEscan_5455MAGs/ISEscan_res/5455MAG_allIS.tsv。注意每个.tsv文件都有标题行，合并时不要重复合并标题行。
# 合并文件（保留第一个文件的标题，跳过后续文件的标题）
awk 'FNR==1 && NR!=1 {next} {print}' \
    ${PROJECT_ROOT:-/path/to/08ESKAPE}/6ISEscan/results/renamed_4538genomes/*.tsv \
    > ${PROJECT_ROOT:-/path/to/08ESKAPE}/6ISEscan/ISEscan_res/4538genomes_allIS.tsv

#只保留5455MAG_allIS.tsv的标题行以及type列为c的行
head -n 1 4538genomes_allIS.tsv | tr "\t" '\n' | nl -v 1 #确认type是第22列

#只保留了完整的IS。序列长度从400bp到3/4000bp，从所有1091IS中筛到309个完整IS。
awk -F"\t" 'NR==1 || $22=="c"' 4538genomes_allIS.tsv > 4538genomes_filtered_completeIS.tsv

#后续分析：提取这些IS所在contig预测DS；删除这些IS后重新预测DS。前后如果存在不同，则说明细菌通过IS调节了DS组成，调节了defensome。

