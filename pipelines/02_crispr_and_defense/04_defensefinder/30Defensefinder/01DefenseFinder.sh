#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J df2
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=112


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 04_defensefinder
# Defense system annotation of the ESKAPE genomes with DefenseFinder
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:30Defensefinder/01DefenseFinder.sh
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

MAGs_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_4538_ESKAPE_genomes_renamed.fa
defense_workdir=${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder
output_dir="${defense_workdir}/res1"

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
