# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 04_spacer_target_genes
# Extraction and annotation of the genes targeted by CRISPR spacers
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:54spacer_target_gene/methylase_phylogeny/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/methylase_phylogeny

tail -n+2 ./ESKAPE_phage_methylase.tsv | cut -f1 | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/merged_spacer_targeted_genes_dereplicated.faa > ./spacer_targeted_methylase.faa

#病毒methylase去冗余
source activate mmseqs2
mmseqs easy-cluster ./spacer_targeted_methylase.faa \
    nr_spacer_targeted_methylase \
	mmseqs2_clustering_protein \
    --min-seq-id 0.6 \
    -c 0.8 \
    --cov-mode 1 \
    --threads 28 \
    -e 1e-5
# 95% 序列一致性
# 85% 覆盖度（AF）
# 覆盖度模式：80% of bothe q and t
# 聚类模式
#CD-HIT applies a greedy incremental clustering strategy (in MMseqs2--cluster-mode 2) to cluster its sequences. MMSeqs2 automatically picks the optimal clustering strategy based on the coverage mode (--cov-mode 0 = set cover,--cov-mode 1,2 = greedy incremental). CD-HIT in default computes the sequence identity globally (-G 1), while MMseqs2 computes it locally. 
# 线程数
# e-value 阈值

#手工整理四个新病毒的methyltransferase蛋白
grep 'vire_genome_803640' ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/VIRE_spacer_targeted_genes_info.txt
echo "vire_genome_803640_30" | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/ESKAPE_phages_from_vire.faa >> ./four_new_phage_methylase.faa

grep 'vire_genome_1727010' ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/VIRE_spacer_targeted_genes_info.txt
echo "vire_genome_1727010_10" | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/merged_spacer_targeted_genes_dereplicated.faa >> ./four_new_phage_methylase.faa

grep 'vire_genome_12877' ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/VIRE_spacer_targeted_genes_info.txt
echo "vire_genome_12877_32" | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/merged_spacer_targeted_genes_dereplicated.faa >> ./four_new_phage_methylase.faa

grep 'IMGVR_UViG_3300035126_000413' ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/IMGVR_spacer_targeted_genes_info.txt
echo "IMGVR_UViG_3300035126_000413|3300035126|Ga0335057_0000171_48" | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/ESKAPE_phages_from_IMGVR.faa >> ./four_new_phage_methylase.faa

cat ./four_new_phage_methylase.faa ./nr_spacer_targeted_methylase_rep_seq.fasta > all_phage_spacer_targeted_methylase.faa

#########提取ESKAPE基因组中的R-M系统中的M############################
grep 'MTase' ${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1/all_4538_ESKAPE_genomes_renamed_defense_finder_genes.tsv > ./host_M.tsv

cut -f2 ./host_M.tsv | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1/all_4538_ESKAPE_genomes_renamed.prt > ./host_M.faa

source activate mmseqs2
mmseqs easy-cluster ./host_M.faa \
    nr_host_M \
	mmseqs2_clustering_protein \
    --min-seq-id 0.6 \
    -c 0.8 \
    --cov-mode 1 \
    --threads 28 \
    -e 1e-5

cat ./all_phage_spacer_targeted_methylase.faa ./nr_host_M_rep_seq.fasta > ./merged_methylase.faa



##提取用于结构分析的蛋白

#type II MTase
seqkit grep -p "GCA_000828935_1-1_2875" ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/methylase_phylogeny/host_M.faa

#type I MTase
seqkit grep -p "GCA_001516105_1-1_1751" ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/methylase_phylogeny/host_M.faa

#type III MTase
seqkit grep -p "GCA_002946935_1-1_721" ${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/methylase_phylogeny/host_M.faa

#直接取序列去AF预测结构
