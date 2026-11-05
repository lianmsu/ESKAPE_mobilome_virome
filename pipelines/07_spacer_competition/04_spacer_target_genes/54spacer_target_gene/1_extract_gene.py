#!/usr/bin/env python3
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 04_spacer_target_genes
# Extraction and annotation of the genes targeted by CRISPR spacers
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:54spacer_target_gene/1_extract_gene.py
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

import re
import sys
from collections import defaultdict

def parse_faa(faa_file):
    """解析faa文件，提取基因位置信息"""
    genes = {} # 创建一个空字典，用于存储所有基因信息
    current_gene = None # 当前正在处理的基因
    current_seq = [] # 当前基因的序列片段列表
    
    with open(faa_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # 保存上一个基因
                if current_gene:
                    genes[current_gene['id']] = {
                        'header': current_gene['header'],
                        'start': current_gene['start'],
                        'end': current_gene['end'],
                        'sequence': ''.join(current_seq)
                    }
                
                # 解析新的基因标题
                # 格式：>CJ_14-10BD_SAB-bin.1-1_1 # 3 # 1406 # 1 # ...
                match = re.match(r'^>([^#]+) # (\d+) # (\d+)', line)
                if match:
                    current_gene = {
                        'id': match.group(1).strip(),
                        'start': int(match.group(2)),
                        'end': int(match.group(3)),
                        'header': line
                    }
                    current_seq = []
                else:
                    current_gene = None
            else:
                if current_gene:
                    current_seq.append(line)
        
        # 保存最后一个基因
        if current_gene:
            genes[current_gene['id']] = {
                'header': current_gene['header'],
                'start': current_gene['start'],
                'end': current_gene['end'],
                'sequence': ''.join(current_seq)
            }
    
    return genes

def parse_protospacers(protospacer_file):
    """解析spacer文件"""
    protospacers = [] # 创建空列表存储所有protospacer信息
    with open(protospacer_file, 'r') as f:
        # 跳过表头
        header = f.readline()
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    protospacers.append({
                        'phage_id': parts[0],
                        'begin': int(parts[1]),
                        'end': int(parts[2])
                    })
    return protospacers

def find_targeted_genes(genes, protospacers):
    """筛选被spacer靶向的基因"""
    targeted = []
    
    for protospacer in protospacers:
        protospacer_begin = protospacer['begin']
        protospacer_end = protospacer['end']
        phage_id = protospacer['phage_id']
        
        # 查找匹配的基因
        for gene_id, gene_info in genes.items():
            # 检查基因ID是否包含phage_id
            if phage_id in gene_id:
                gene_start = gene_info['start']
                gene_end = gene_info['end']
                
                # 检查protospacer是否在基因范围内
                if gene_start <= min(protospacer_begin, protospacer_end) and gene_end >= min(protospacer_begin, protospacer_end):
                    targeted.append({
                        'gene_id': gene_id,
                        'gene_start': gene_start,
                        'gene_end': gene_end,
                        'protospacer_phage': phage_id,
                        'protospacer_begin': protospacer_begin,
                        'protospacer_end': protospacer_end,
                        'gene_info': gene_info
                    })
                    break
    
    return targeted

def main():
    # 文件路径
    protospacer_file = "${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/spacers_IMGVR.blastn.1mis100cov_phage_filtered.tsv"
    faa_file = "${PROJECT_ROOT:-/path/to/08ESKAPE}/54spacer_target_gene/ESKAPE_phages_from_IMGVR.faa"  # 请替换为实际的faa文件路径
    output_file = "IMGVR_spacer_targeted_genes.faa"
    info_file = "IMGVR_spacer_targeted_genes_info.txt"
    
    print("解析faa文件...")
    genes = parse_faa(faa_file)
    print(f"找到 {len(genes)} 个基因")
    
    print("解析spacer文件...")
    protospacers = parse_protospacers(protospacer_file)
    print(f"找到 {len(protospacers)} 个protospacer")
    
    print("筛选靶向基因...")
    targeted = find_targeted_genes(genes, protospacers)
    print(f"找到 {len(targeted)} 个被靶向的基因")
    
    # 保存被靶向的基因序列
    with open(output_file, 'w') as f:
        for item in targeted:
            f.write(item['gene_info']['header'] + '\n')
            # 每60个字符换行
            seq = item['gene_info']['sequence']
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')
    
    # 保存详细信息
    with open(info_file, 'w') as f:
        f.write("Gene_ID\tGene_Start\tGene_End\tProtospacer_Phage\tProtospacer_Begin\tProtospacer_End\n")
        for item in targeted:
            f.write(f"{item['gene_id']}\t{item['gene_start']}\t{item['gene_end']}\t"
                   f"{item['protospacer_phage']}\t{item['protospacer_begin']}\t{item['protospacer_end']}\n")
    
    print(f"\n结果已保存:")
    print(f"  - 基因序列: {output_file}")
    print(f"  - 详细信息: {info_file}")

if __name__ == "__main__":
    main()
