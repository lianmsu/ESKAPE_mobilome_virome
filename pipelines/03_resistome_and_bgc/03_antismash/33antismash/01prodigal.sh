#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J prodigal
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.01prodigal_lpw.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 03_resistome_and_bgc / 03_antismash
# Biosynthetic gene cluster prediction with antiSMASH
# Manuscript      : Fig. 1D
# Original path   : lfat1t01:33antismash/01prodigal.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate prokka
# 开始计时
start_time=$(date +%s)

cd ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/41antiSMASH

prodigal -p meta -i ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/30minced_5455MAG/merged_5455MAG.fa -f gbk -o 5455MAG.gbk -a 5455MAGproteins.faa -d 5455MAGgenes.fna -g 11 -q


# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
