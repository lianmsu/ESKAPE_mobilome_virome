#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J integron_4538
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=450GB


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 01_ice_integron_isescan
# Integron detection with IntegronFinder
# Manuscript      : Fig. 2A
# Original path   : pkusz_CPU:4integronfinder/01integron_finder.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

# 开始计时
start_time=$(date +%s)



#########################################################################
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/4integronfinder
source activate integron_finder
#########################################################################


mags="${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep/all_4538_ESKAPE_genomes_renamed.fa"

integron_finder --local-max --func-annot ${mags} --cpu $SLURM_NTASKS --outdir ${PROJECT_ROOT:-/path/to/08ESKAPE}/4integronfinder/res --gbk --pdf


# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
