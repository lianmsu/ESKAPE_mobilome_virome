
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 05_virus_discovery / 04_final_catalogue
# Filtering rules that define the final non-redundant ESKAPE virome
# Manuscript      : Fig. 3
# Original path   : lfat1t01:46final_virus/1virus_filter.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus
#对spacer匹配到的病毒进行筛选，长度>5kb，For geNomad predictions, contigs of 5-10 kb were required to have a virus score ≥ 0.9, at least one viral hallmark gene, and a virus-marker enrichment > 2.0, whereas contigs ≥ 10 kb were required to have a virus score ≥ 0.8 and at least one hallmark gene or a virus-marker enrichment > 5.0. Other contigs with direct or inverted terminal repeats (DTRs or ITRs) were also maintained.

head -n1 ${PROJECT_ROOT:-/path/to/08ESKAPE}/45virus_genomad/res4/checkv_clustered_modified_phage_summary/checkv_clustered_modified_phage_virus_summary.tsv > genomad_id_MG.tsv
awk -F "\t" 'NR>1 {
    # 5-10kb: score>=0.9, hallmarks>=1, enrichment>2.0
    if (($3 ~ /DTR|ITR/) || 
		($2 >= 5000 && $2 < 10000 && $7 >= 0.9 && $9 >= 1 && $10 > 2.0) || 
    # >=10kb: score>=0.8, (hallmarks>=1 OR enrichment>5.0)
        ($2 >= 10000 && $7 >= 0.8 && ($9 >= 1 || $10 > 5.0))) 
    print $0
}' ${PROJECT_ROOT:-/path/to/08ESKAPE}/45virus_genomad/res4/checkv_clustered_modified_phage_summary/checkv_clustered_modified_phage_virus_summary.tsv >> genomad_id_MG.tsv

# 匹配 |provirus_数字_数字 结尾的模式
# 去除 |provirus 及其后的所有内容（包括数字等）
cut -f1 ./genomad_id_MG.tsv | sed 's/|provirus.*$//' > ./genomad_clean_id_MG.txt

#
#checkv viral_gene>0
head -n1 ${PROJECT_ROOT:-/path/to/08ESKAPE}/45virus_checkv/quality_summary.tsv > checkv_virus_gene1.tsv
awk -F "\t" 'NR>1 {
	if($6>=1)
	print $0
}' ${PROJECT_ROOT:-/path/to/08ESKAPE}/45virus_checkv/quality_summary.tsv >> checkv_virus_gene1.tsv

#
# 一行命令完成
awk 'NR==FNR{a[$1]; next} $1 in a' genomad_clean_id_MG.txt checkv_virus_gene1.tsv | cut -f1 > filtered_phage.id.txt

seqkit grep -f filtered_phage.id.txt ${PROJECT_ROOT:-/path/to/08ESKAPE}/44virus_checkv_cluster/checkv_clustered_modified_phage.fa > ./checkv_clustered_genomad_checkv_filtered_phage.fa
#16905

#从checkv_virus_gene1.tsv和genomad_id_MG.tsv中提取对应的详细质量和物种信息
#checkv_virus_gene1.tsv的序列名和输入序列一致，软件没有加额外后缀。
#genomad_id_MG.tsv的序列名可能加了|provirus后缀。
#所以对genomad_id_MG.tsv的第一列删除|provirus后缀，注意这个文件的长度信息不能用
#只保留两个文件中的标题行以及序列名包含filtered_phage.id.txt中序列名的行。
#!/bin/bash

# 1. 清理genomad文件的ID
awk 'BEGIN{FS=OFS="\t"} NR==1{print; next} {sub(/\|provirus.*/, "", $1); print}' genomad_id_MG.tsv > genomad_clean.tsv

# 2. 提取两个文件的对应行
# 先获取两个文件的标题
checkv_header=$(head -1 checkv_virus_gene1.tsv)
genomad_header=$(head -1 genomad_clean.tsv)

# 创建输出文件标题
echo -e "ID\t${checkv_header}\t${genomad_header}" | sed "s/\t\t/\t/g" > phage_detailed_info.tsv

# 3. 对每个ID进行提取
while read id; do
    # 从checkv文件提取
    checkv_line=$(grep -w "^$id" checkv_virus_gene1.tsv)
    
    # 从genomad文件提取
    genomad_line=$(grep -w "^$id" genomad_clean.tsv)
    
    if [ -n "$checkv_line" ] && [ -n "$genomad_line" ]; then
        # 合并两行，用ID开头
        echo -e "${id}\t${checkv_line}\t${genomad_line}" >> phage_detailed_info.tsv
   # elif [ -n "$checkv_line" ]; then
        # 只在checkv中找到
   #     genomad_empty=$(printf '\t%.0s' {1..$(awk 'NR==1{print NF-1}' genomad_clean.tsv)})
   #     echo -e "${id}\t${checkv_line}${genomad_empty}" >> phage_detailed_info.tsv
   # elif [ -n "$genomad_line" ]; then
        # 只在genomad中找到
   #     checkv_empty=$(printf '\t%.0s' {1..$(awk 'NR==1{print NF-1}' checkv_virus_gene1.tsv)})
   #     echo -e "${id}\t${checkv_empty}${genomad_line}" >> phage_detailed_info.tsv
    fi
done < filtered_phage.id.txt

echo "完成！结果保存在 phage_detailed_info.tsv"

head -1 phage_detailed_info.tsv > phage_detailed_info_trim.tsv
awk -F "\t" 'NR>1 && NF == 26' phage_detailed_info.tsv >> phage_detailed_info_trim.tsv


