#!/usr/bin/env python3

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 09_phylogeny / 01_marker_genes
# hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL
# Manuscript      : Fig. 6
# Original path   : lfat1t01:60marker_gene_phylogeny/2_anno_treatment_32phage.py
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

import os
import glob
import re

# 设置输入和输出目录
input_dir = "32phage_hhsuite"
output_dir = "${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/32phage_anno"

# 创建输出目录
os.makedirs(output_dir, exist_ok=True)

# 获取所有.hhr文件
hhr_files = glob.glob(os.path.join(input_dir, "*.hhr"))

for hhr_file in hhr_files:
    # 获取基本文件名
    base_name = os.path.basename(hhr_file).replace('.hhr', '')
    tsv_file = os.path.join(input_dir, f"{base_name}.tsv")
    output_file = os.path.join(output_dir, f"{base_name}_filtered.hhr")
    
    print(f"处理文件: {base_name}")
    
    # 检查.tsv文件是否存在
    if not os.path.exists(tsv_file):
        print(f"警告: 找不到对应的.tsv文件: {tsv_file}")
        continue
    
    # 读取.tsv文件内容
    with open(tsv_file, 'r') as f:
        tsv_lines = [line.strip().replace('>', '') for line in f.readlines() if line.strip()]
    
    # 读取.hhr文件内容
    with open(hhr_file, 'r') as f:
        hhr_content = f.read()
    
    # 找到表格开始和结束的位置
    table_start = hhr_content.find("No Hit")
    table_end = hhr_content.find("No 1")
    
    if table_start == -1:
        print(f"错误: 在文件 {hhr_file} 中找不到表格开始位置")
        continue
    
    if table_end == -1:
        # 如果没有"No 1"，则取到文件末尾
        table_end = len(hhr_content)
    
    # 提取表格内容
    table_content = hhr_content[table_start:table_end].strip()
    table_lines = table_content.split('\n')
    
    # 写入输出文件
    with open(output_file, 'w') as out_f:
        # 写入表头行（添加detail_anno列）
        header_line = table_lines[0].strip()
        out_f.write(header_line + "  detail_anno\n")
        
        # 处理数据行
        tsv_index = 0
        for line in table_lines[1:]:
            line = line.strip()
            if not line:
                continue
            
            # 检查是否是数据行（以数字开头）
            if re.match(r'^\s*\d+', line):
                # 使用更精确的方法提取Prob值
                # 假设Prob值在行的固定位置（通常在Hit名称后的第一个数字）
                parts = line.split()
                
                # 找到第一个看起来像概率的数字（0-100之间的数字）
                prob = None
                for i, part in enumerate(parts):
                    if i < 3:  # 跳过前几个字段（序号、Hit名称）
                        continue
                    try:
                        if '.' in part and 0 <= float(part) <= 100:
                            prob = float(part)
                            break
                    except ValueError:
                        continue
                
                if prob is not None and prob > 90 and tsv_index < len(tsv_lines):
                    # 获取对应的.tsv注释
                    detail_anno = tsv_lines[tsv_index]
                    # 写入合并后的行
                    out_f.write(f"{line}  {detail_anno}\n")
                
                tsv_index += 1
    
    print(f"完成: {output_file}")

print("所有文件处理完成！")
