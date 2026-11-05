#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J df2
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 04_virus_defense
# Defense system annotation of the viruses
# Manuscript      : Fig. 3F
# Original path   : lfat1t01:47virus_defensefinder/01DefenseFinder.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

# 开始计时
start_time=$(date +%s)

source activate defensefinder
#########################################################################

#########################################################################

MAGs_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.fa
defense_workdir=${PROJECT_ROOT:-/path/to/08ESKAPE}/47virus_defensefinder
output_dir="${defense_workdir}/res2"

mkdir -p ${output_dir}
cd ${defense_workdir}

defense-finder run -a ${MAGs_path} --preserve-raw -o "$output_dir"


# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
