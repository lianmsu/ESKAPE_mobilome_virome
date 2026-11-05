#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J prodigal
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=112
#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 06_virus_and_host_databases
# VIRE download and ORF prediction for the ESKAPE virus references
# Manuscript      : Fig. 4B
# Original path   : lfat1t01:2virus_database/2VIRE/1prodigal.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate instrain
# 开始计时
start_time=$(date +%s)

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE
mags=${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE/vire_all_genomes_combined.fasta
output_dir=${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE

echo "Running Prodigal"
prodigal -p meta -i "${mags}" \
		-f gff -o "${output_dir}/vire_all_genomes_combined.gff" \
		-m -d "${output_dir}/vire_all_genomes_combined.fna" \
		-a "${output_dir}/vire_all_genomes_combined.faa"

grep -v '^#' vire_all_genomes_combined.gff > vire_all_genomes_combined_filtered.gff

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
