# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 10_experiments / 01_phage_sequence
# Codon usage and sequence analysis of the isolated phages
# Manuscript      : Fig. 6
# Original path   : lfat1t01:70experiments/01phage_sequence/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/01phage_sequence

seqkit grep -p "vire_genome_12877" ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/potential_novel_virus.fa > ./vire_genome_12877.fa

seqkit grep -p "vire_genome_1727010" ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/potential_novel_virus.fa > ./vire_genome_1727010.fa

seqkit grep -p "vire_genome_803640" ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/potential_novel_virus.fa > ./vire_genome_803640.fa

seqkit grep -p "IMGVR_UViG_3300035126_000413|3300035126|Ga0335057_0000171" ${PROJECT_ROOT:-/path/to/08ESKAPE}/52potential_novel_virus/potential_novel_virus.fa > ./IMGVR_UViG_3300035126_000413.fa


#去除边缘的repeat，然后用prodigal的nromal模式预测基因
mamba activate python3.12.4
python3 apc.py -b vire_genome_12877 vire_genome_12877.fa
python3 apc.py -b vire_genome_1727010 vire_genome_1727010.fa
python3 apc.py -b vire_genome_803640 vire_genome_803640.fa
python3 apc.py -b IMGVR_UViG_3300035126_000413 IMGVR_UViG_3300035126_000413.fa



mamba activate prokka
cat ./*trimmed* > all_4_novel_phage.fa


prodigal -p single -i ./all_4_novel_phage.fa -f gff -o all_4_novel_phage.gff -a all_4_novel_phage.faa -d all_4_novel_phage.fna -g 11 -n -q

#删除.gff文件中以#开头的行
# 单个文件处理
#grep -v '^#' vire_genome_12877_trimmed.1.gff > vire_genome_12877_trimmed.1.gff

