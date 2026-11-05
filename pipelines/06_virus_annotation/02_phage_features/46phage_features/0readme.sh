# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 02_phage_features
# Molecular features of the phages: GC, NARSC, length, quality
# Manuscript      : Fig. 3
# Original path   : lfat1t01:46phage_features/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_features

#统计病毒的基因组特性，GC 长度 蛋白性质等

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_features

mkdir -p genomes
mkdir -p proteins

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.fa ./

# 输出：序列ID \t 长度 \t GC含量 \t 序列(可选)
seqkit fx2tab -n -l -g -H ${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_features/ESKAPE_all_phage_info_filter.fa > per_phage_stats.tsv

#complete phages
grep 'Complete' ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.tsv | cut -f1 | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/46iPHoP/ESKAPE_all_phage_info_filter.fa > ESKAPE_complete_phage.fa

seqkit split -i ESKAPE_complete_phage.fa -O ./genomes/

cd genomes
for f in ESKAPE_complete_phage.part_*.fa; do
    # 提取 "part_" 之后的所有内容作为新名字
    newname="${f#*part_}"
    
    # 只有当新名字和旧名字不同时才执行移动，避免报错
    if [ "$f" != "$newname" ]; then
        mv "$f" "$newname"
    fi
done

source activate prokka

# 设置路径
GENOME_DIR="${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_features/genomes"
PROTEIN_DIR="${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_features/proteins"  # 或者指定其他路径

# 创建输出目录
mkdir -p "$PROTEIN_DIR"

# 计数器
count=0
total=$(ls "$GENOME_DIR"/*.fa 2>/dev/null | wc -l)

echo "开始处理 $total 个基因组文件..."

# 遍历所有 .fa 文件
for genome_file in "$GENOME_DIR"/*.fa; do
    # 获取文件名（不带路径）
    filename=$(basename "$genome_file")
    
    # 生成输出文件名（将 .fa 替换为 .faa）
    output_file="${filename%.fa}.faa"
    output_path="$PROTEIN_DIR/$output_file"
    
    # 运行 Prodigal
    echo "处理第 $((++count))/$total 个文件: $filename"
    
    # 使用 -a 参数输出蛋白质序列，-p meta 模式适用于宏基因组数据
    prodigal -i "$genome_file" -a "$output_path" -p meta -q
    
    # 或者使用更详细的参数，根据你的数据类型选择：
    # prodigal -i "$genome_file" -a "$output_path" -d "${PROTEIN_DIR}/${filename%.fa}.ffn" -o "${PROTEIN_DIR}/${filename%.fa}.gff" -p meta -q
    
    # 检查 Prodigal 是否成功
    if [ $? -eq 0 ]; then
        echo "  √ 成功: $output_file"
    else
        echo "  × 失败: $filename"
    fi
done

