#!/usr/bin/env python3

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 09_phylogeny / 01_marker_genes
# hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL
# Manuscript      : Fig. 6
# Original path   : lfat1t01:60marker_gene_phylogeny/3_anno_summary_32phage.py
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

import os
import glob
import re
import pandas as pd

# 设置目录路径
input_dir = "${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/32phage_anno"
output_file = "32phage_annotation_summary.tsv"

# 获取所有_filtered.hhr文件
file_pattern = os.path.join(input_dir, "*_filtered.hhr")
files = glob.glob(file_pattern)

# 存储结果的列表
results = []

for file_path in files:
    # 提取文件名
    file_name = os.path.basename(file_path).replace('_filtered.hhr', '')
    
    # 初始化默认值
    prob = "unannotated"
    detail_anno = "unannotated"
    
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            
            # 检查文件是否有数据行
            if len(lines) >= 2:
                data_line = lines[1].strip()
                
                if data_line:
                    # 使用正则表达式查找Prob值（格式如 99.7）
                    prob_match = re.search(r'(\d+\.\d+)\s+[\d\.E\-\+]+', data_line)
                    
                    if prob_match:
                        prob = prob_match.group(1)
                        # 提取detail_anno部分（Prob之后的所有内容）
                        detail_start = prob_match.end()
                        detail_anno = data_line[detail_start:].strip()
                        
                        # 如果detail_anno以数字或科学计数法开头，可能是E-value等，需要进一步处理
                        # 找到第一个非数字的字段作为detail_anno的开始
                        detail_parts = detail_anno.split()
                        for i, part in enumerate(detail_parts):
                            if not re.match(r'^[\d\.E\-\+]+$', part):
                                detail_anno = ' '.join(detail_parts[i:])
                                break
    
    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {e}")
        prob = "error"
        detail_anno = f"error: {str(e)}"
    
    results.append({
        'seq': file_name,
        'Prob': prob,
        'detail_anno': detail_anno
    })

# 创建DataFrame并保存
df = pd.DataFrame(results)
df.to_csv(output_file, sep='\t', index=False)

print(f"处理完成！共处理 {len(results)} 个文件。")
print(f"结果已保存到: {output_file}")
print("\n前10行结果预览:")
print(df.head(10))
