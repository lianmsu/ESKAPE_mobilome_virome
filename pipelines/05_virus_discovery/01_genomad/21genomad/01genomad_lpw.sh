#!/bin/sh
#SBATCH --partition=cpu
#SBATCH -J gmMAG
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH -o job.genomad_lpw.%j.%N.out

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 05_virus_discovery / 01_genomad
# geNomad virus and plasmid detection on the ESKAPE genomes
# Manuscript      : Fig. 3
# Original path   : lfat1t01:21genomad/01genomad_lpw.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

# 开始计时
start_time=$(date +%s)


source activate genomad
genomad_db=${DB_ROOT:-/path/to/databases}/genomad_db
#########################################################################
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/21genomad
MAG_contig=${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_4538_ESKAPE_genomes_renamed.fa
workdir=${PROJECT_ROOT:-/path/to/08ESKAPE}/21genomad

genomad end-to-end --cleanup ${MAG_contig} ${workdir}/res ${genomad_db}

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "脚本运行时间：$runtime 秒"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0



