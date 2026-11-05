
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 10_experiments / 04_strain_prophage
# Prophage content of the experimental host strains
# Original path   : lfat1t01:70experiments/04strain_prophage/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/04strain_prophage

#实验所用菌株中是否有prophage，DS


#Kp菌株ATCC 13883
GCF_000742135.1
GCF_000788015.1

mamba activate ncbi-genome-download
datasets download genome accession --inputfile kp_ATCC13883_accession --dehydrated --filename kp.zip
unzip kp.zip -d kp
datasets rehydrate --directory kp/

mkdir -p ./all_genomes
find ${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/04strain_prophage/kp/ncbi_dataset -type f -name "*.fna" -exec sh -c '
  for f do
    base=$(basename "$f")          # 去除前缀

    cp "$f" "./all_genomes/$base"
  done
' sh {} +


#Ab菌株ATCC 19606
GCF_009035845.1
GCF_020911985.1
GCF_014672755.1
GCF_009759685.1
GCF_025995075.1
GCF_019331655.1

datasets download genome accession --inputfile Ab_ATCC19606_accession --dehydrated --filename ab.zip
unzip ab.zip -d ab
datasets rehydrate --directory ab/

find ${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/04strain_prophage/ab/ncbi_dataset -type f -name "*.fna" -exec sh -c '
  for f do
    base=$(basename "$f")          # 去除前缀

    cp "$f" "./all_genomes/$base"
  done
' sh {} +

cd ./all_genomes
for file in *.fna; do
    newname=$(echo "$file" | awk -F'_' '{print $1"_"$2".fa"}')
    #echo "测试模式: $file → $newname"
    # 实际执行时取消下一行的注释
    mv -v "$file" "$newname"
done

for file in *.fa; do
    [ -e "$file" ] || continue
    
    # 分离文件名主体和后缀
    # basename "$file" .fa 会去掉 .fa 后缀，只留下主体部分
    name_body=$(basename "$file" .fa)
    
    # 只替换主体部分中的 . 为 _
    new_body="${name_body//./_}"
    
    # 重组新文件名
    newname="${new_body}.fa"
    
    if [ "$file" != "$newname" ]; then
        #echo "重命名: $file -> $newname"
        mv -v "$file" "$newname"
    fi
done

for file in *.fa; do
    prefix=$(basename "$file" .fa)
    tmpfile="${file}.tmp"
    awk -v prefix="$prefix" '
    BEGIN {counter=1}
    /^>/ {
        if($0 ~ /plasmid/) {
            printf ">%s-%d-plasmid\n", prefix, counter++;
        } else {
            printf ">%s-%d\n", prefix, counter++;
        }
        next
    }
    {print}' "$file" > "$tmpfile"
    mv "$tmpfile" "$file"
done

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/04strain_prophage
cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/04strain_prophage/all_genomes/*.fa > merged_genome.fa

mamba activate genomad
genomad_db=${DB_ROOT:-/path/to/databases}/genomad_db
genomad end-to-end --cleanup merged_genome.fa ./genomad_res ${genomad_db}

source activate checkv
export CHECKVDB=${DB_ROOT:-/path/to/databases}/checkv_db/checkv-db-v1.5
virus_contig=${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/04strain_prophage/genomad_res/merged_genome_summary/merged_genome_virus.fna
workdir=${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/04strain_prophage/checkv
#end-to-end：污染度、完整度分析和去除原噬菌体污染
checkv end_to_end ${virus_contig} ${workdir} -t 24

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/04strain_prophage
source activate defensefinder
defense-finder run -a ./merged_genome.fa --preserve-raw -o ./DF_res

####################################################################
mkdir -p att_site

#目的：找到attP序列，看待验证的温和噬菌体上有没有一致的attP序列。如果有，可能会被抑制。
