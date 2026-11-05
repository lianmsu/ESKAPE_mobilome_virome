#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J is4538
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=450GB
#SBATCH --output=job.lpw-isescan4538.%j.%N.out


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 01_ice_integron_isescan
# Insertion sequence detection with ISEScan
# Manuscript      : Fig. 2A
# Original path   : pkusz_CPU:6ISEscan/ISEScan_4538MAGs.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate isescan
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/6ISEscan
id_file=${PROJECT_ROOT:-/path/to/08ESKAPE}/6ISEscan/4538MAG_path.txt
mkdir -p ./results
#start=15001
#end=20000
# 切片
#ids=$(awk "NR >= $start && NR <= $end" ${id_file})

#for i in ${ids}
#do

#isescan.py --seqfile ${HOME_ROOT:-/path/to/home}/03GW/01ANI99MAGs/${i} --output ./results/${i} --nthread 24

#done

cat 4538MAG_path.txt | xargs -n 1 -P 8 -I{} isescan.py --seqfile {} --output ./results --nthread 6

