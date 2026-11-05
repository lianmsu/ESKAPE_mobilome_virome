# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 00_notes
# Defense system annotation of the ESKAPE genomes with DefenseFinder
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:30Defensefinder/res1/DG_gff/readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1/DG_gff

#根据文件${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/res1/merged_5467MAG_defense_finder_genes.tsv的第二列从${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/12prodigal_5467MAGs/merged_5467MAG_ORF.gff中提取对应的行到一个新文件DG_gff.gff

# 定义输入和输出文件路径
DEFENSE_FINDER_FILE="${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1/all_4538_ESKAPE_genomes_renamed_defense_finder_genes.tsv"
GFF_FILE="${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/12prodigal_5467MAGs/merged_5467MAG_ORF.gff"
OUTPUT_FILE="DG_gff.gff"

# # 检查输入文件是否存在
# if [ ! -f "$DEFENSE_FINDER_FILE" ] || [ ! -f "$GFF_FILE" ]; then
    # echo "错误：输入文件不存在"
    # exit 1
# fi

# # 使用awk处理文件
# awk -F"\t" '
# BEGIN {
    # # 初始化计数器
    # processed_lines = 0
    # matched_lines = 0
# }

# # 处理DefenseFinder文件（第一个文件）
# NR==FNR {
    # # 跳过标题行
    # if (NR == 1) next
    
    # # 存储hit_id到数组，使用哈希表提高查找效率
    # if ($2 != "") {
        # hit_ids[$2] = 1
    # }
    # next
# }

# # 处理GFF文件（第二个文件）
# {
    # # 跳过注释行和空行
    # if ($0 ~ /^#/ || $0 ~ /^$/) next
    
    # processed_lines++
    
    # # 从第一列获取序列ID
    # seq_id = $1
    
    # # 从第九列提取ID值
    # if (match($9, /ID=([^;]+)/, m)) {
        # # 提取基因ID
        # gene_id = m[1]
        
        # # 分割基因ID并获取最后一部分
        # split(gene_id, arr, "_")
        # id = arr[length(arr)]
        
        # # 构建组合ID
        # combined_id = seq_id "_" id
        # echo $combined_id
        # # 检查是否存在于hit_ids中
        # if (combined_id in hit_ids) {
            # print $0
            # matched_lines++
        # }
    # }
# }

# END {
    # # 输出处理统计信息
    # printf "处理完成：\n"
    # printf "总处理行数：%d\n", processed_lines
    # printf "匹配行数：%d\n", matched_lines
# }
# ' "$DEFENSE_FINDER_FILE" "$GFF_FILE" > "$OUTPUT_FILE"

# # 检查输出文件是否成功创建
# if [ -f "$OUTPUT_FILE" ]; then
    # echo "处理完成，结果已保存到 $OUTPUT_FILE"
# else
    # echo "错误：输出文件创建失败"
    # exit 1
# fi 






# 定义输入和输出文件路径
DEFENSE_FINDER_FILE="${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1/all_4538_ESKAPE_genomes_renamed_defense_finder_genes.tsv"
GFF_FILE="${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1/all_4538_ESKAPE_genomes_renamed.prt.gff"
OUTPUT_FILE="DF_gene_gff.gff"

# 检查输入文件是否存在
if [ ! -f "$DEFENSE_FINDER_FILE" ] || [ ! -f "$GFF_FILE" ]; then
    echo "错误：输入文件不存在"
    exit 1
fi

# 使用awk处理文件
awk -F"\t" '
BEGIN {
    # 初始化计数器
    processed_lines = 0
    matched_lines = 0
}

# 处理DefenseFinder文件（第一个文件）
NR==FNR {
    # 跳过标题行
    if (NR == 1) next
    
    # 存储hit_id到数组，使用哈希表提高查找效率
    if ($2 != "") {
        hit_ids[$2] = 1
    }
    next
}

# 处理GFF文件（第二个文件）
{
    # 跳过注释行和空行
    if ($0 ~ /^#/ || $0 ~ /^$/) next
    
    processed_lines++
    
    # 从第一列获取序列ID
    seq_id = $1
    
    # 从第九列提取ID值
    if (match($9, /ID=([^;]+)/, m)) {
        # 提取基因ID
        gene_id = m[1]
        
        # 分割基因ID并获取最后一部分
        split(gene_id, arr, "_")
        id = arr[length(arr)]
        
        # 构建组合ID
        combined_id = seq_id "_" id
        echo $combined_id
        # 检查是否存在于hit_ids中
        if (combined_id in hit_ids) {
            print $0
            matched_lines++
        }
    }
}

END {
    # 输出处理统计信息
    printf "处理完成：\n"
    printf "总处理行数：%d\n", processed_lines
    printf "匹配行数：%d\n", matched_lines
}
' "$DEFENSE_FINDER_FILE" "$GFF_FILE" > "$OUTPUT_FILE"

# 检查输出文件是否成功创建
if [ -f "$OUTPUT_FILE" ]; then
    echo "处理完成，结果已保存到 $OUTPUT_FILE"
else
    echo "错误：输出文件创建失败"
    exit 1
fi
