#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J cv
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=end


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 05_virus_discovery / 02_checkv
# CheckV on the second virus batch
# Manuscript      : Fig. 3A
# Original path   : lfat1t01:45virus_checkv/01checkv_lpw.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

#进入checkv并导入数据库目录
source activate checkv
export CHECKVDB=${DB_ROOT:-/path/to/databases}/checkv_db/checkv-db-v1.5

virus_contig=${PROJECT_ROOT:-/path/to/08ESKAPE}/44virus_checkv_cluster/checkv_clustered_modified_phage.fa
workdir=${PROJECT_ROOT:-/path/to/08ESKAPE}/45virus_checkv
#end-to-end：污染度、完整度分析和去除原噬菌体污染
checkv end_to_end ${virus_contig} ${workdir} -t $SLURM_NTASKS
