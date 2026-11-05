# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 09_phylogeny / 01_marker_genes
# hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL
# Manuscript      : Fig. 6
# Original path   : lfat1t01:60marker_gene_phylogeny/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny

#所有893个完整ESKAPE病毒的MCP terL protal序列系统发育##############################

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/50vcontact3/ESKAPE_complete_phage.fa ./

mamba activate prokka
prodigal -p meta -i ESKAPE_complete_phage.fa \
		-f gff -o ./ESKAPE_complete_phage.gff \
		-m -d ./ESKAPE_complete_phage.fna \
		-a ./ESKAPE_complete_phage.faa

grep -v '^#' ESKAPE_complete_phage.gff > ESKAPE_complete_phage_treated.gff

#hhsuite注释
#本地整理ORF与注释的关系
#挑选目标基因进行建树


#额外32个病毒序列的标志基因也要加到分析中#############################################
cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/32_phage.fasta ./

mamba activate prokka
prodigal -p meta -i 32_phage.fasta \
		-f gff -o ./32_phage.gff \
		-m -d ./32_phage.fna \
		-a ./32_phage.faa

grep -v '^#' 32_phage.gff > 32_phage_treated.gff
