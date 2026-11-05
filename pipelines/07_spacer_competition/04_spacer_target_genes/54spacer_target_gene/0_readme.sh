# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 04_spacer_target_genes
# Extraction and annotation of the genes targeted by CRISPR spacers
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:54spacer_target_gene/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene

#提取被spacer靶向的病毒基因，对这些基因进行注释，看主要是哪些基因

#1. 来自IMGVR的病毒
ln -s ${PROJECT_ROOT:-/path/to/08ESKAPE}/40spacer2IMGVR/spacers_IMGVR.blastn.1mis100cov.out ./
#-outfmt '6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send'

cut -f2,15,16 ./spacers_IMGVR.blastn.1mis100cov.out | sort | uniq > spacers_IMGVR.blastn.1mis100cov_phage.tsv

awk 'BEGIN {FS=OFS="\t"} 
     NR==FNR {ids[$1]; next} 
     ($1 in ids)' \
     ESKAPE_all_phage_info_filter.tsv \
     spacers_IMGVR.blastn.1mis100cov_phage.tsv \
     > spacers_IMGVR.blastn.1mis100cov_phage_filtered.tsv
# 直接修改原文件（添加表头）
sed -i "1i phageID\tprotospacer_begin\tprotospacer_end" spacers_IMGVR.blastn.1mis100cov_phage_filtered.tsv

cut -f1 ESKAPE_all_phage_info_filter.tsv | grep 'IMGVR' | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.fa > ./ESKAPE_phages_from_IMGVR.fa

source activate prokka
prodigal -p meta -i ./ESKAPE_phages_from_IMGVR.fa \
    -f gff -o ./ESKAPE_phages_from_IMGVR.gff \
    -m -d ./ESKAPE_phages_from_IMGVR.fna \
    -a ./ESKAPE_phages_from_IMGVR.faa

#根据protospacer的位置信息和gff文件中的位置信息来筛选被spacer靶向的基因
#运行python代码完成

########################################################################
#1. 来自VIRE的病毒
ln -s ${PROJECT_ROOT:-/path/to/08ESKAPE}/41spacer2VIRE/spacers_VIRE.blastn.1mis100cov.out ./
#-outfmt '6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send'

cut -f2,15,16 ./spacers_VIRE.blastn.1mis100cov.out | sort | uniq > spacers_VIRE.blastn.1mis100cov_phage.tsv

awk 'BEGIN {FS=OFS="\t"} 
     NR==FNR {ids[$1]; next} 
     ($1 in ids)' \
     ESKAPE_all_phage_info_filter.tsv \
     spacers_VIRE.blastn.1mis100cov_phage.tsv \
     > spacers_VIRE.blastn.1mis100cov_phage_filtered.tsv
# 直接修改原文件（添加表头）
sed -i "1i phageID\tprotospacer_begin\tprotospacer_end" spacers_VIRE.blastn.1mis100cov_phage_filtered.tsv

cut -f1 ESKAPE_all_phage_info_filter.tsv | grep 'vire' | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.fa > ./ESKAPE_phages_from_vire.fa

source activate prokka
prodigal -p meta -i ./ESKAPE_phages_from_vire.fa \
    -f gff -o ./ESKAPE_phages_from_vire.gff \
    -m -d ./ESKAPE_phages_from_vire.fna \
    -a ./ESKAPE_phages_from_vire.faa

#根据protospacer的位置信息和gff文件中的位置信息来筛选被spacer靶向的基因
#运行python代码完成

#除了这些含有protospacer的基因，还有一些基因含有多个protospacer。

tail -n+2 ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/IMGVR_spacer_targeted_genes_info.txt | cut -f1 | sort | uniq | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/ESKAPE_phages_from_IMGVR.faa > IMGVR_spacer_targeted_genes_dereplicated.faa

tail -n+2 ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/VIRE_spacer_targeted_genes_info.txt | cut -f1 | sort | uniq | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/ESKAPE_phages_from_vire.faa > VIRE_spacer_targeted_genes_dereplicated.faa

#对得到的基因进行注释 eggnog
cat IMGVR_spacer_targeted_genes_dereplicated.faa VIRE_spacer_targeted_genes_dereplicated.faa > merged_spacer_targeted_genes_dereplicated.faa
