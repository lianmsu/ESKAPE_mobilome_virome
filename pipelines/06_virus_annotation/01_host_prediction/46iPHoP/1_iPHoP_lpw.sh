#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J iphop_esk
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 01_host_prediction
# Host prediction with iPHoP and lifestyle with BACPHLIP
# Manuscript      : Fig. 3C-E
# Original path   : lfat1t01:46iPHoP/1_iPHoP_lpw.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

# 记录开始时间（秒级时间戳）
start_time=$(date +%s)

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP
my_input_phages=${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus/checkv_clustered_genomad_checkv_filtered_phage.fa
db_path=${DB_ROOT:-/path/to/databases}/iPHoP_db/Jun_2025_pub_rw/
out_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/iphop_output_refdb

source activate iphopv142

iphop predict --fa_file ${my_input_phages} --db_dir ${db_path} --out_dir ${out_path} -t 56 --min_score 90 --add_tsv_output

end_time=$(date +%s)
cost_time=$(( end_time - start_time ))
echo "execution time is $(($cost_time/60))min $(($cost_time%60))s"
