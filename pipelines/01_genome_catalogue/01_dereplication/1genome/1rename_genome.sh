# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 01_genome_catalogue / 01_dereplication
# Download and rename the complete ESKAPE genomes from NCBI RefSeq
# Manuscript      : Fig. 1A
# Original path   : pkusz_CPU:1genome/1rename_genome.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome

cd ./all_genomes

for file in *.fna; do
    newname=$(echo "$file" | awk -F'_' '{print $1"_"$2".fa"}')
    #echo "测试模式: $file → $newname"
    # 实际执行时取消下一行的注释
    mv -v "$file" "$newname"
done

#contig重命名
input_dir="${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_genomes"

for input_file in "$input_dir"/*.fa; do
    basename=$(basename "$input_file" .fa)
    awk -v basename="$basename" '/^>/ {print ">" basename "-" ++i; next} {print}' "$input_file" > tmpfile && mv tmpfile "$input_file"
    echo "已处理: $input_file"
done


cat ./*.fa > all_xxx_ESKAPE_genomes

sed '/^>/ s/\./_/g' all_4696_predators_genome.fa > all_4696_predators_genome_renamed.fa
