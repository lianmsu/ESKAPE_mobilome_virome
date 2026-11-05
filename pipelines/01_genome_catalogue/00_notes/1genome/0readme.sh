
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 01_genome_catalogue / 00_notes
# Genome renaming, GTDB metadata matching and plasmid extraction commands
# Manuscript      : Fig. 1
# Original path   : lfat1t01:1genome/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/renamed_4538genomes
for file in *.fa; do
sed -i '/^>/ s/\./_/g' ${file}

done

#得到基因组的其他信息，根据文件ESKAPE_genome_info.tsv的第一列，在文件${HOME_ROOT:-/path/to/home}/06mobilome_predators/1genome_collections/1GTDB_R226/bac120_metadata_r226.tsv中进行匹配，提取匹配到的行
# 提取ID列表（跳过可能的表头）
tail -n +2 ESKAPE_genome_info.tsv | cut -f1 > genome_ids.txt

# 使用grep匹配整行（-w确保精确匹配整个字段）
grep -w -f genome_ids.txt ${HOME_ROOT:-/path/to/home}/06mobilome_predators/1genome_collections/1GTDB_R226/bac120_metadata_r226.tsv > genome_matched_GTDB_metadata.tsv

# 提取ID列表（跳过可能的表头）
tail -n +2 ESKAPE_genome_info.tsv | cut -f17 > biosample_ids.txt

# 使用grep匹配整行（-w确保精确匹配整个字段）
grep -w -f biosample_ids.txt ${HOME_ROOT:-/path/to/home}/06mobilome_predators/1genome_collections/1GTDB_R226/bac120_metadata_r226.tsv > biosample_matched_GTDB_metadata.tsv

############提取质粒信息############
seqkit fx2tab \
  --name \
  --length \
  --gc \
  --header-line \
  ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_4538_ESKAPE_genomes_renamed.fa | \
  grep -i "plasmid" > ./plasmid_info.tsv


