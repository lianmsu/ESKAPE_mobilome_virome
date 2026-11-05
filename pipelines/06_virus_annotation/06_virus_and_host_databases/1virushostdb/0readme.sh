# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 06_virus_and_host_databases
# ESKAPE virus subset of Virus-Host DB
# Manuscript      : Fig. 4A
# Original path   : lfat1t01:2virus_database/1virushostdb/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb

gunzip -d virushostdb.genomic.fna.gz

grep 'Enterococcus faecium' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.tsv | cut -f4 | wc -l
#51

grep 'Staphylococcus aureus' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.tsv | cut -f4 | wc -l
#512

grep 'Klebsiella pneumoniae' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.tsv | cut -f4 | wc -l
#1475

grep 'Acinetobacter baumannii' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.tsv | cut -f4 | wc -l
#340

grep 'Pseudomonas aeruginosa' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.tsv | cut -f4 | wc -l
#1066

grep 'Enterobacter ' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.tsv | cut -f4 | wc -l
#121

seqkit grep -n -r -i -p "Enterococcus faecium" ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.formatted.genomic.fna -o Ef_vhdb.fasta
#50
seqkit grep -n -r -i -p "Staphylococcus aureus" ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.formatted.genomic.fna -o Sa_vhdb.fasta
#535
seqkit grep -n -r -i -p "Klebsiella pneumoniae" ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.formatted.genomic.fna -o Kp_vhdb.fasta
#1490
seqkit grep -n -r -i -p "Acinetobacter baumannii" ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.formatted.genomic.fna -o Ab_vhdb.fasta
#345
seqkit grep -n -r -i -p "Pseudomonas aeruginosa" ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.formatted.genomic.fna -o Pa_vhdb.fasta
#1046
seqkit grep -n -r -i -p "Enterobacter " ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.formatted.genomic.fna -o En_vhdb.fasta
#124

cat ./*.fasta > vhdb_ESKAPE_virus.fasta
# 直接在原文件上修改（小心！建议先备份）
sed -i 's/^>\([^ ]*\).*/>\1/' vhdb_ESKAPE_virus.fasta
#3590个病毒序列

#先提取FASTA序列名，再用这些名字去TSV文件中匹配。
# 第一步：从FASTA提取所有序列名（去掉'>'）
grep '^>' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/vhdb_ESKAPE_virus.fasta | sed 's/^>//' > sequence_names.txt

# 第二步：使用grep -f匹配TSV文件中的这些名字
head -n1 ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.tsv > VHDB_ESKAPE_virus.tsv
grep -f sequence_names.txt ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/1virushostdb/virushostdb.tsv >> VHDB_ESKAPE_virus.tsv

