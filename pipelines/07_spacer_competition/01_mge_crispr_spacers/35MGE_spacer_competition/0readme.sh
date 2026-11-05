# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 01_mge_crispr_spacers
# Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets)
# Manuscript      : Fig. 2B
# Original path   : lfat1t01:35MGE_spacer_competition/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition

#整理各种MGE对应序列，分别注释其中的spacer，把spacer掩蔽后进行blastn

#1、提取plasmid
grep 'plasmid' ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_4538_ESKAPE_genomes_renamed.fa | \
  seqkit seq -n | \
  seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_4538_ESKAPE_genomes_renamed.fa > ./all_plasmid_contig.fa

#2、提取provirus#
cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/42genomad_checkv/ESKAPE_genome_prophage_filtered.fa ./

#3、integron序列
# 创建BED格式的区域文件
tail -n +2 ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/ESKAPE_integron_n.tsv | \
  awk '{print $2"\t"$3"\t"$4}' > integron_regions.bed
# 使用BED文件提取序列
seqkit subseq --bed integron_regions.bed \
  ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_4538_ESKAPE_genomes_renamed.fa \
  -o ESKAPE_complete_integrons.fa
# 清理临时文件（可选）
rm integron_regions.bed

# 直接修改FASTA文件
sed "/^>/ s/:\.\s*$//" ./ESKAPE_complete_integrons.fa \
  > ESKAPE_complete_integrons_clean.fa

#ICE 序列
#从深圳服务器复制



