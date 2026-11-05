# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 08_phage_clustering_and_host / 03_novel_virus
# Identification and comparison of the potential novel viruses
# Manuscript      : Fig. 5-6
# Original path   : lfat1t01:52potential_novel_virus/clinker_virus/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus

grep 'vire_genome_803640' ${PROJECT_ROOT:-/path/to/08ESKAPE}/44virus_checkv_cluster/my_clusters.tsv 

grep 'IMGVR_UViG_3300035126_000413' ${PROJECT_ROOT:-/path/to/08ESKAPE}/44virus_checkv_cluster/my_clusters.tsv 
IMGVR_UViG_3300035126_000413|3300035126|Ga0335057_0000171	IMGVR_UViG_3300035126_000413|3300035126|Ga0335057_0000171,IMGVR_UViG_3300044729_001474|3300044729|Ga0453742_0000480,IMGVR_UViG_3300044638_003266|3300044638|Ga0453716_0000317,vire_genome_149183,IMGVR_UViG_3300043464_001031|3300043464|Ga0453100_0000520,IMGVR_UViG_3300042528_000127|3300042528|Ga0453971_0000203,IMGVR_UViG_3300034062_003382|3300034062|Ga0334995_0005768,IMGVR_UViG_3300044675_005077|3300044675|Ga0453754_0006433,IMGVR_UViG_3300034092_002356|3300034092|Ga0335010_0012705,IMGVR_UViG_3300034111_001199|3300034111|Ga0335063_0014533

grep 'vire_genome_1727010' ${PROJECT_ROOT:-/path/to/08ESKAPE}/44virus_checkv_cluster/my_clusters.tsv 

grep 'vire_genome_12877' ${PROJECT_ROOT:-/path/to/08ESKAPE}/44virus_checkv_cluster/my_clusters.tsv 

##因为寻找新病毒时，vcontact3只用了完整的病毒基因组，
#要用它们的近源不完整序列做genomic comparison

#用所有病毒序列进行vcontact3
grep 'vire_genome_803640' ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus/vcontact3_res_all_phage/exports/networks/part1_cosmo_data.csv

grep 'IMGVR_UViG_3300035126_000413' ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus/vcontact3_res_all_phage/exports/networks/part1_cosmo_data.csv

grep 'vire_genome_1727010' ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus/vcontact3_res_all_phage/exports/networks/part1_cosmo_data.csv

grep 'vire_genome_12877' ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus/vcontact3_res_all_phage/exports/networks/part1_cosmo_data.csv



#手动整理

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus
echo 'IMGVR_UViG_3300035126_000413|3300035126|Ga0335057_0000171,IMGVR_UViG_3300044729_001474|3300044729|Ga0453742_0000480,IMGVR_UViG_3300044638_003266|3300044638|Ga0453716_0000317,vire_genome_149183,IMGVR_UViG_3300043464_001031|3300043464|Ga0453100_0000520,IMGVR_UViG_3300042528_000127|3300042528|Ga0453971_0000203,IMGVR_UViG_3300034062_003382|3300034062|Ga0334995_0005768,IMGVR_UViG_3300044675_005077|3300044675|Ga0453754_0006433,IMGVR_UViG_3300034092_002356|3300034092|Ga0335010_0012705,IMGVR_UViG_3300034111_001199|3300034111|Ga0335063_0014533' | tr ',' '\n' | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/43combined_all_virus/all_ESKAPE_phage.fa > Ab_v_2_groups.fa

#
mkdir -p ./phage_sequence

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus/selected_phage_for_clinker.txt | while read -r i ; do

seqkit grep -p "${i}" ${PROJECT_ROOT:-/path/to/08ESKAPE}/43combined_all_virus/all_ESKAPE_phage.fa > ./phage_sequence/${i}.fasta

done


#得到这些质粒的hhsuite注释
#查找含有'IMGPR_plasmid_3300009181_000028|3300009181|Ga0114969_10015712'的蛋白序列在文件${HOME_ROOT:-/path/to/home}/06mobilome_predators/80HHsuit/nr_plasmid_rep_seq.faa中的顺序i，根据顺序i提取${HOME_ROOT:-/path/to/home}/06mobilome_predators/80HHsuit/nr_plasmid_hhsuite下的seq{i}.*到路径${HOME_ROOT:-/path/to/home}/06mobilome_predators/80HHsuit/special_plasmid/IMGPR_plasmid_3300009181_000028


#用antismash为三个质粒基因组生成gbk文件

mamba activate antiSMASH

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus/selected_phage_for_clinker.txt | while read -r i ; do

antismash ./phage_sequence/${i}.fasta \
	--taxon bacteria \
	--cpus 10 \
	--databases ${HOME_ROOT:-/path/to/home}/software/miniforge3/envs/antiSMASH/lib/python3.11/site-packages/antismash/databases \
	--output-dir ./${i} \
	--allow-long-headers \
	--genefinding-tool prodigal-m

done

####################################
mkdir -p ./phage_sequence

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus/selected_phage_for_clinker2.txt | while read -r i ; do

seqkit grep -p "${i}" ${PROJECT_ROOT:-/path/to/08ESKAPE}/43combined_all_virus/all_ESKAPE_phage.fa > ./phage_sequence/${i}.fasta

done

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/clinker_virus/selected_phage_for_clinker2.txt | while read -r i ; do

antismash ./phage_sequence/${i}.fasta \
	--taxon bacteria \
	--cpus 10 \
	--databases ${HOME_ROOT:-/path/to/home}/software/miniforge3/envs/antiSMASH/lib/python3.11/site-packages/antismash/databases \
	--output-dir ./${i} \
	--allow-long-headers \
	--genefinding-tool prodigal-m

done
