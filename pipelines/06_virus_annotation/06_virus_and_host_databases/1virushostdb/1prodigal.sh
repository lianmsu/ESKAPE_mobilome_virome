#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J prodigal
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 06_virus_and_host_databases
# ESKAPE virus subset of Virus-Host DB
# Manuscript      : Fig. 4A
# Original path   : lfat1t01:2virus_database/1virushostdb/1prodigal.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate instrain
# 开始计时
start_time=$(date +%s)

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb
mags=${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.genomic.fna
output_dir=${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb

echo "Running Prodigal"
prodigal -p meta -i "${mags}" \
		-f gff -o "${output_dir}/virushostdb.gff" \
		-m -d "${output_dir}/virushostdb.fna" \
		-a "${output_dir}/virushostdb.faa"

grep -v '^#' virushostdb.gff > virushostdb_filtered.gff

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
