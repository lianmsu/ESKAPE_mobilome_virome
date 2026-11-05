# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 08_phage_clustering_and_host / 04_multi_host_phage
# Multi-host phages and their host ranges
# Manuscript      : Fig. 4
# Original path   : lfat1t01:53multi_host_phage/clinker/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/53multi_host_phage/clinker


#手动整理

echo 'vire_genome_39529,IMGVR_UViG_3300048922_000045|3300048922|Ga0496119_0000589,vire_genome_871096,vire_genome_1377468,vire_genome_93879,vire_genome_444250' | tr ',' '\n' | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.fa > selected_6phages.fa

#
mkdir -p ./phage_sequence

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/53multi_host_phage/clinker/selected_phage_for_clinker.txt | while read -r i ; do

seqkit grep -p "${i}" ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.fa > ./phage_sequence/${i}.fasta

done


#得到这些质粒的hhsuite注释
#查找含有'IMGPR_plasmid_3300009181_000028|3300009181|Ga0114969_10015712'的蛋白序列在文件${HOME_ROOT:-/path/to/home}/06mobilome_predators/80HHsuit/nr_plasmid_rep_seq.faa中的顺序i，根据顺序i提取${HOME_ROOT:-/path/to/home}/06mobilome_predators/80HHsuit/nr_plasmid_hhsuite下的seq{i}.*到路径${HOME_ROOT:-/path/to/home}/06mobilome_predators/80HHsuit/special_plasmid/IMGPR_plasmid_3300009181_000028


#用antismash为三个质粒基因组生成gbk文件

mamba activate antiSMASH

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/53multi_host_phage/clinker/selected_phage_for_clinker.txt | while read -r i ; do

antismash ./phage_sequence/${i}.fasta \
	--taxon bacteria \
	--cpus 10 \
	--databases ${HOME_ROOT:-/path/to/home}/software/miniforge3/envs/antiSMASH/lib/python3.11/site-packages/antismash/databases \
	--output-dir ./${i} \
	--allow-long-headers \
	--genefinding-tool prodigal-m

done

####################################
source activate prokka
prodigal -p meta -i ./selected_6phages.fa \
    -f gff -o ./selected_6phages.gff \
    -m -d ./selected_6phages.fna \
    -a ./selected_6phages.faa
