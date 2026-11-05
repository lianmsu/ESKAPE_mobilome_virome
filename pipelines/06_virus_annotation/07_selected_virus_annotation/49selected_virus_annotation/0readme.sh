# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 07_selected_virus_annotation
# hhsuite2 annotation of the selected (novel) viruses
# Manuscript      : Fig. 5-6
# Original path   : lfat1t01:49selected_virus_annotation/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/49selected_virus_annotation

cut -f1 ${PROJECT_ROOT:-/path/to/08ESKAPE}/49selected_virus_annotation/selected_phage_ID.tsv | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus/ESKAPE_complete_phage.fa > ./selected_phage.fa

source activate prokka
prodigal -p meta -i ./selected_phage.fa -f gff -o selected_phage.gff -a selected_phage.faa -d selected_phage.fna -g 11 -n -q

#对蛋白进行hhsuite注释


#删除.gff文件中以#开头的行
# 单个文件处理
grep -v '^#' selected_phage.gff > selected_phage_treated.gff

# 在本地用R语言处理文件Candidate_17_new_phages_treated.gff，得到proksee需要的格式

# 统计hhsuite注释结果，方便人工检查和确定


seqkit fx2tab -l -g -n -i -H Candidate_17_new_phages.fasta > Candidate_17_new_phages_GC.tsv




