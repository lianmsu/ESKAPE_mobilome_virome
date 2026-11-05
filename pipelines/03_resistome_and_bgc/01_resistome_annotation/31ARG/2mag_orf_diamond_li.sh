#!/bin/sh
#SBATCH -o jobdia1060.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J diamond
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 03_resistome_and_bgc / 01_resistome_annotation
# ARG annotation: RGI/CARD and diamond against SARG
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:31ARG/2mag_orf_diamond_li.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate eggnog

orf="${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG/4538genomes_prodigal.faa"
out=${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG

#各种类型功能基因的数据库
ARGdb=${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/42ARG/database/ARG/CARD.dmnd
MRGdb=${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/42ARG/database/MRG/MRG.dmnd
POPDegdb=${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/42ARG/database/mibPOP/mibPOP.dmnd
VFGdb=${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/42ARG/database/VFG/VFDB_core.dmnd

query=${orf}
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG
mkdir -p ./1ARG
mkdir -p ./2VFG
mkdir -p ./3MRG
mkdir -p ./4mibPOP

diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${ARGdb} -q ${query} -o ${out}/1ARG/diamond_ARG_out --id 60 --subject-cover 80 --query-cover 80
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${MRGdb} -q ${query} -o ${out}/3MRG/diamond_MRG_out --id 60 --subject-cover 80 --query-cover 80
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${POPDegdb} -q ${query} -o ${out}/4mibPOP/diamond_POPDeg_out --id 60 --subject-cover 80 --query-cover 80
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${VFGdb} -q ${query} -o ${out}/2VFG/diamond_VFG_out --id 60 --subject-cover 80 --query-cover 80
#--id 80 比较好
#The host system is detected to have 404 GB of RAM. It is recommended to increase the block size for better performance using these parameters : -b8 -c1
#--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
#

#--id 80, 
echo "successfully"

