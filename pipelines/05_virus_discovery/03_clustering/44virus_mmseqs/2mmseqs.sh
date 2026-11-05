#!/bin/sh
#SBATCH -o job-mmseq2-%j.%N.out
#SBATCH --partition=fat
#SBATCH -J mmseq2
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 05_virus_discovery / 03_clustering
# mmseqs2 clustering of the virus set
# Manuscript      : Fig. 3
# Original path   : lfat1t01:44virus_mmseqs/2mmseqs.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate mmseqs2
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/44virus_mmseqs
# 创建输出目录

mkdir -p mmseqs2_clustering_phage

# 使用 easy-cluster 进行聚类
mmseqs easy-cluster ../43combined_all_virus/all_ESKAPE_phage.fa nr_phage mmseqs2_clustering_phage \
    --min-seq-id 0.95 \
    -c 0.9 \
    --cov-mode 0 \
    --cluster-mode 1 \
    --threads 48 \
    -e 1e-5
# 95% 序列一致性
# 90% 覆盖度（AF）
# 覆盖度模式：90% of bothe q and t
# 贪婪聚类
# 线程数
# e-value 阈值

