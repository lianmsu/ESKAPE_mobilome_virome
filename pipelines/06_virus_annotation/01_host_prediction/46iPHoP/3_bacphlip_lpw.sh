#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J bacph
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 01_host_prediction
# Host prediction with iPHoP and lifestyle with BACPHLIP
# Manuscript      : Fig. 3C-E
# Original path   : lfat1t01:46iPHoP/3_bacphlip_lpw.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

# 记录开始时间（秒级时间戳）
start_time=$(date +%s)

#进入checkv并导入数据库目录
source activate BACPHLIP

virus_contig=${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.fa
workdir=${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP

cd ${workdir}
bacphlip -i ${virus_contig} --multi_fasta


end_time=$(date +%s)
cost_time=$(( end_time - start_time ))
echo "execution time is $(($cost_time/60))min $(($cost_time%60))s"


