#!/bin/sh
#SBATCH -o job.tree.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J MCP
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 09_phylogeny / 01_marker_genes
# hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL
# Manuscript      : Fig. 6
# Original path   : lfat1t01:60marker_gene_phylogeny/MCP/1_iqtree_39phageMCP.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/MCP
source activate iqtree2

#1、先用muscle把protein序列对齐形成多序列比对（MSA）
#muscle -align all_MCP.faa -output all_MCP_mafft.afa -threads 48
linsi --thread 25 ./39phage_MCP.faa > 39phage_MCP_mafft.afa
#mafft --auto --thread 24 383acr.faa > 383acr_mafft.afa

#2、然后用trimal修剪
#trimal -in all_MCP_mafft.afa -gt 0.5 > all_MCP_aln.trim.afa
#-gt 0.5代表删去有50%以上gap的位点。-gt 0.5这个是对应这篇文章方法中的参数 https://www.nature.com/articles/s41467-023-38400-0

#clipkit ${special_gene}_mafft.afa -m kpic-gappy -g 0.9
#-gt 0.5代表删去有50%以上gap的位点。-gt 0.5这个是对应这篇文章方法中的参数 https://www.nature.com/articles/s41467-023-38400-0

clipkit 39phage_MCP_mafft.afa -m gappy -g 0.9


#3、对修建后的MSA进行建树
iqtree2 -s 39phage_MCP_mafft.afa.clipkit -m MFP+MERGE -B 1000 -nt AUTO --nmax 4000

echo "${special_gene} tree finised!"

