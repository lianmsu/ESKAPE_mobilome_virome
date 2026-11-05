# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 03_crisprcaster
# CRISPR-Cas subtyping with CRISPRCasTyper
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:24CRISPRCasTyper/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/24CRISPRCasTyper

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/23CRISPRCasFinder/CRISPRCasFinder-master/Result_minced_crispr_contig_1766546889/TSV/Crisprs_REPORT.tsv ./

sed -i '/^$/d' ./Crisprs_REPORT.tsv

#level4的CRISPR存在于同一contig上的情况
awk -F"\t" '{
    count[$2]++; 
    line[$2]=line[$2] ? line[$2] ORS $0 : $0
} 
END {
    for (id in count)
        if (count[id]>1) print line[id]
}' ./Crisprs_REPORT.tsv


# 1. 提取ID并去重
cut -f2 ./Crisprs_REPORT.tsv | 
sort -u > crispr_ids.tmp
# 2. 筛选序列（使用8线程）
seqkit grep -j 8 -f crispr_ids.tmp \
    ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_4538_ESKAPE_genomes_renamed.fa \
    > ./all_CRISPR_level4.fa

# 3. 验证结果
seqkit stats ./all_CRISPR_level4.fa

#只保留含有Cas operon的contig进行下一步的NR过滤。

