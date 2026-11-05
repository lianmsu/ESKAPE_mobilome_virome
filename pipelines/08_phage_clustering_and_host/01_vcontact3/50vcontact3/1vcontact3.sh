#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J vc3
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 08_phage_clustering_and_host / 01_vcontact3
# vContact3 clustering of the ESKAPE viruses
# Manuscript      : Fig. 5
# Original path   : lfat1t01:50vcontact3/1vcontact3.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

# 开始计时
start_time=$(date +%s)

source activate vcontact3
#########################################################################
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/50vcontact3
#########################################################################
phage_genome_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/46final_virus/ESKAPE_all_phage.fa

vcontact3 run --nucleotide ${phage_genome_path} \
	--output results_dir2 \
	--db-domain "prokaryotes" \
	--db-path ${DB_ROOT:-/path/to/databases}/vcontact3/latest \
	--db-version 230 \
	--exports {cytoscape,graphml,cosmograph,d3js,ani,newick,profiles,completeness,centroids} \
	--breaks 10 \
	--threads 56

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
