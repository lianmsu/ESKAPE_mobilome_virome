#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J eggnog
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 04_spacer_target_genes
# Extraction and annotation of the genes targeted by CRISPR spacers
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:54spacer_target_gene/2_eggnog_mapper.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate eggnog

input_proteins=${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/merged_spacer_targeted_genes_dereplicated.faa
output_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/eggnog_res
database_path=${DB_ROOT:-/path/to/databases}/eggnog/

mkdir -p ${output_path}
cd ${output_path}

# 开始计时
start_time=$(date +%s)
#echo "当前目录为：$(pwd)"

#参数来自eggnog-mapper网页版的默认参数
emapper.py --cpu $SLURM_NTASKS --data_dir ${database_path} -o out --output_dir ${output_path} --temp_dir ${output_path} --override -m diamond --dmnd_ignore_warnings -i ${input_proteins} --evalue 0.00001 --score 60 --pident 60 --query_cover 80 --subject_cover 80 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel
##删除了--dmnd_algo ctg 

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
