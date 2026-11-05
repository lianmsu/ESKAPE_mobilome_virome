#!/usr/bin/env python3
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 01_ice_integron_isescan
# ICE detection with ICEfinder2 (batches of 1000 and the full set)
# Manuscript      : Fig. 2A
# Original path   : pkusz_CPU:5ICEfinder/ICEres/process_iceid.py
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

import re

def process_ice_sequences(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    results = []
    for line in lines:
        # 移除'>'符号
        line = line.strip().lstrip('>')
        
        # 使用正则表达式提取contig id
        contig_match = re.search(r'^.*?contig_\d+', line)
        if not contig_match:
            continue
        contig_id = contig_match.group()
        
        # 提取位置信息
        position_match = re.search(r'(\d+)-(\d+)$', line)
        if not position_match:
            continue
        start, end = map(int, position_match.groups())
        
        # 计算序列长度
        length = end - start + 1
        
        # 格式化输出
        result = f"{contig_id}\t{start}-{end}\t{length}"
        results.append(result)
    
    # 写入结果到输出文件
    with open(output_file, 'w') as f:
        for result in results:
            f.write(result + '\n')

if __name__ == "__main__":
    input_file = "all_ICE_contigid.txt"
    output_file = "processed_ice_sequences_ID.txt"
    process_ice_sequences(input_file, output_file)
