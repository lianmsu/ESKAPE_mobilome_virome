#!/bin/sh
#SBATCH -o job-ICE_4538-%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J ICE_4538
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=250GB


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 01_ice_integron_isescan
# ICE detection with ICEfinder2 (batches of 1000 and the full set)
# Manuscript      : Fig. 2A
# Original path   : pkusz_CPU:5ICEfinder/icefinder_4538.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

# 开始计时
start_time=$(date +%s)
source activate icefinder
#########################################################################
# File Name: GO_DefenseFinder.sh
# Date: 17/4/2024
#########################################################################

cp -r ${HOME_ROOT:-/path/to/home}/05plateau/21icefinder_5455MAGs/ICEfinder2_moban ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEfinder2_4538
ICE_workdir=${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/ICEfinder2_4538
cd ${ICE_workdir}
id_file=${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder/4538MAG.id
MAG_dir=${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep/renamed_4538genomes

start=3001
end=4538
# 切片
ids=$(awk "NR >= $start && NR <= $end" ${id_file})

for i in ${ids}
do

base=$(basename "$i" .fa)
python ICEfinder2.py -i ${MAG_dir}/"${base}_renamed.fa" -t Metagenome

done

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
