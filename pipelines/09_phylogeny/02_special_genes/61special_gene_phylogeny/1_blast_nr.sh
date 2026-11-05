#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J blast
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.blast.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 09_phylogeny / 02_special_genes
# Phylogenies of the special genes (e.g. VanX) found in the new viruses
# Manuscript      : Fig. 6
# Original path   : lfat1t01:61special_gene_phylogeny/1_blast_nr.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/61special_gene_phylogeny
# 设置BLAST数据库路径
source activate blast2.16
export BLASTDB=${DB_ROOT:-/path/to/databases}/NR/nr
# 蛋白序列使用blastp
blastp -query ./phage_Dps.faa \
       -db nr \
	   -task blastp \
       -out Dps_vs_nr_blast.raw \
       -outfmt "6 std staxids stitle" \
       -evalue 1e-5 \
       -mt_mode 1 \
       -num_threads 28 \
       -max_target_seqs 100

blastp -query ./phage_VanX.faa \
       -db nr \
       -task blastp \
       -out VanX_vs_nr_blast.raw \
       -outfmt "6 std staxids stitle" \
       -evalue 1e-5 \
       -mt_mode 1 \
       -num_threads 28 \
       -max_target_seqs 100
