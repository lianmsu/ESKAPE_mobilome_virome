# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 02_phage_features
# Molecular features of the phages: GC, NARSC, length, quality
# Manuscript      : Fig. 3
# Original path   : lfat1t01:46phage_features/get_gc_and_narsc_python3.py
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

# -*- coding: utf-8 -*-
"""
Updated for Python 3.x
Comprehensive Molecular Characteristics Scanner
Includes: GC, Nc, GC_rank, NARSC, CARSC, SARSC, AvgMW
"""

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# 自动处理 Biopython GC 计算函数迁移
try:
    from Bio.SeqUtils import gc_fraction
except ImportError:
    from Bio.SeqUtils import GC as gc_fraction

# ==================== 1. 数据字典补全 ====================

# 包含 U, J, B, Z 的完整氨基酸原子字典
aa_dictionary = {
    'K': {'N': 1, 'S': 0, 'C': 4, 'MW': 146.1882}, 'R': {'N': 3, 'S': 0, 'C': 4, 'MW': 174.2017},
    'H': {'N': 2, 'S': 0, 'C': 4, 'MW': 155.1552}, 'D': {'N': 0, 'S': 0, 'C': 2, 'MW': 133.1032},
    'E': {'N': 0, 'S': 0, 'C': 3, 'MW': 147.1299}, 'N': {'N': 1, 'S': 0, 'C': 2, 'MW': 132.1184},
    'Q': {'N': 1, 'S': 0, 'C': 3, 'MW': 146.1451}, 'S': {'N': 0, 'S': 0, 'C': 1, 'MW': 105.0930},
    'T': {'N': 0, 'S': 0, 'C': 2, 'MW': 119.1197}, 'Y': {'N': 0, 'S': 0, 'C': 7, 'MW': 181.1894},
    'A': {'N': 0, 'S': 0, 'C': 1, 'MW': 89.0935},  'V': {'N': 0, 'S': 0, 'C': 3, 'MW': 117.1469},
    'L': {'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736}, 'I': {'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736},
    'P': {'N': 0, 'S': 0, 'C': 3, 'MW': 115.1310}, 'F': {'N': 0, 'S': 0, 'C': 7, 'MW': 165.1900},
    'M': {'N': 0, 'S': 1, 'C': 3, 'MW': 149.2124}, 'W': {'N': 1, 'S': 0, 'C': 9, 'MW': 204.2262},
    'G': {'N': 0, 'S': 0, 'C': 0, 'MW': 75.0669},  'C': {'N': 0, 'S': 1, 'C': 1, 'MW': 121.1590},
    # 补全特殊编码
    'U': {'N': 0, 'S': 0, 'C': 1, 'MW': 168.07},   # 硒半胱氨酸 (侧链硒取代硫)
    'J': {'N': 0, 'S': 0, 'C': 4, 'MW': 131.1736}, # I/L 简并
    'B': {'N': 0.5, 'S': 0, 'C': 2, 'MW': 132.6108}, # D/N 简并，氮取均值
    'Z': {'N': 0.5, 'S': 0, 'C': 3, 'MW': 146.6375}  # E/Q 简并，氮取均值
}

# 密码子字典（用于计算 Nc 和 GC_rank）
codon_dictionary = {
    'GAT': {'aa': 'D', 'GC_rank': 1}, 'GAC': {'aa': 'D', 'GC_rank': 0},
    'GAA': {'aa': 'E', 'GC_rank': 1}, 'GAG': {'aa': 'E', 'GC_rank': 0},
    'TCT': {'aa': 'S', 'GC_rank': 1}, 'TCC': {'aa': 'S', 'GC_rank': 0},
    'TCA': {'aa': 'S', 'GC_rank': 1}, 'TCG': {'aa': 'S', 'GC_rank': 0},
    'AGT': {'aa': 'S', 'GC_rank': 1}, 'AGC': {'aa': 'S', 'GC_rank': 0},
    'ACT': {'aa': 'T', 'GC_rank': 1}, 'ACC': {'aa': 'T', 'GC_rank': 0},
    'ACA': {'aa': 'T', 'GC_rank': 1}, 'ACG': {'aa': 'T', 'GC_rank': 0},
    'TAT': {'aa': 'Y', 'GC_rank': 1}, 'TAC': {'aa': 'Y', 'GC_rank': 0},
    'GCT': {'aa': 'A', 'GC_rank': 1}, 'GCC': {'aa': 'A', 'GC_rank': 0},
    'GCA': {'aa': 'A', 'GC_rank': 1}, 'GCG': {'aa': 'A', 'GC_rank': 0},
    'GTT': {'aa': 'V', 'GC_rank': 0}, 'GTC': {'aa': 'V', 'GC_rank': 1},
    'GTA': {'aa': 'V', 'GC_rank': 0}, 'GTG': {'aa': 'V', 'GC_rank': 1},
    'TTA': {'aa': 'L', 'GC_rank': 1}, 'TTG': {'aa': 'L', 'GC_rank': .5},
    'CTT': {'aa': 'L', 'GC_rank': .5}, 'CTC': {'aa': 'L', 'GC_rank': 0},
    'CTA': {'aa': 'L', 'GC_rank': .5}, 'CTG': {'aa': 'L', 'GC_rank': 0},
    'ATT': {'aa': 'I', 'GC_rank': 1}, 'ATC': {'aa': 'I', 'GC_rank': 0},
    'ATA': {'aa': 'I', 'GC_rank': 1}, 'CCT': {'aa': 'P', 'GC_rank': 1},
    'CCC': {'aa': 'P', 'GC_rank': 0}, 'CCA': {'aa': 'P', 'GC_rank': 1},
    'CCG': {'aa': 'P', 'GC_rank': 0}, 'TTT': {'aa': 'F', 'GC_rank': 1},
    'TTC': {'aa': 'F', 'GC_rank': 0}, 'GGT': {'aa': 'G', 'GC_rank': 1},
    'GGC': {'aa': 'G', 'GC_rank': 0}, 'GGA': {'aa': 'G', 'GC_rank': 1},
    'GGG': {'aa': 'G', 'GC_rank': 0}, 'TGT': {'aa': 'C', 'GC_rank': 1},
    'TGC': {'aa': 'C', 'GC_rank': 0}, 'ATG': {'aa': 'M', 'GC_rank': None},
    'AAA': {'aa': 'K', 'GC_rank': 1}, 'AAG': {'aa': 'K', 'GC_rank': 0},
    'TGG': {'aa': 'W', 'GC_rank': None}, 'AAT': {'aa': 'N', 'GC_rank': 1},
    'AAC': {'aa': 'N', 'GC_rank': 0}, 'CAA': {'aa': 'Q', 'GC_rank': 1},
    'CAG': {'aa': 'Q', 'GC_rank': 0}, 'CAT': {'aa': 'H', 'GC_rank': 1},
    'CAC': {'aa': 'H', 'GC_rank': 0}, 'CGT': {'aa': 'R', 'GC_rank': .5},
    'CGC': {'aa': 'R', 'GC_rank': 0}, 'CGA': {'aa': 'R', 'GC_rank': .5},
    'CGG': {'aa': 'R', 'GC_rank': 0}, 'AGA': {'aa': 'R', 'GC_rank': 1},
    'AGG': {'aa': 'R', 'GC_rank': .5}
}

# ==================== 2. 计算函数 ====================

def calculate_NC_and_Rank(dna_seq):
    """计算 Nc 和 GC_rank"""
    code_11_aas = {
        'M': ['ATG'], 'W': ['TGG'], 'F': ['TTT', 'TTC'], 'Y': ['TAT', 'TAC'],
        'H': ['CAT', 'CAC'], 'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'],
        'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'],
        'C': ['TGT', 'TGC'], 'I': ['ATT', 'ATC', 'ATA'], 'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'L': ['TTA', 'TTG', 'CTT','CTC','CTA','CTG'],
        'S': ['TCT', 'TCC', 'TCA','TCG','AGT','AGC'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG','AGA','AGG']
    }
    sf_map = {'M':'SF1','W':'SF1','F':'SF2','Y':'SF2','H':'SF2','Q':'SF2','N':'SF2','K':'SF2','D':'SF2','E':'SF2','C':'SF2','I':'SF3','A':'SF4','P':'SF4','T':'SF4','V':'SF4','G':'SF4','L':'SF6','S':'SF6','R':'SF6'}

    counts = {c: 0 for c in codon_dictionary.keys()}
    aa_counts = {aa: 0 for aa in code_11_aas.keys()}
    sum_rank, rank_pos = 0.0, 0
    seq = str(dna_seq).upper()

    for i in range(0, len(seq)-len(seq)%3, 3):
        codon = seq[i:i+3]
        if codon in codon_dictionary:
            counts[codon] += 1
            aa = codon_dictionary[codon]['aa']
            aa_counts[aa] += 1
            if codon_dictionary[codon]['GC_rank'] is not None:
                sum_rank += codon_dictionary[codon]['GC_rank']
                rank_pos += 1

    # 计算 Ne
    sf_groups = {f'SF{i}': [] for i in [1, 2, 3, 4, 6]}
    for aa, codons in code_11_aas.items():
        n = float(aa_counts[aa])
        if n > 1:
            sum_p_sq = sum((counts[c]/n)**2 for c in codons)
            f_val = (n * sum_p_sq - 1) / (n - 1)
            if f_val > 0:
                sf_groups[sf_map[aa]].append(1.0/f_val)

    nc = 2.0
    for sf, counts_list in sf_groups.items():
        if counts_list:
            if sf == 'SF3': nc += counts_list[0]
            elif sf in ['SF2', 'SF4', 'SF6']:
                av_f = sum(1.0/ne for ne in counts_list) / len(counts_list)
                nc += (9 if sf=='SF2' else 5 if sf=='SF4' else 3) / av_f

    avg_rank = round(sum_rank/rank_pos, 4) if rank_pos > 0 else "Nan"
    return round(min(nc, 61.0), 4), avg_rank

def calculate_protein_metrics(aa_seq):
    """计算 NARSC, CARSC, SARSC, AvgMW"""
    clean = [x for x in str(aa_seq).upper() if x in aa_dictionary]
    if not clean: return "0", "0", "0", "0"
    L = len(clean)
    n = sum(aa_dictionary[x]['N'] for x in clean) / L
    c = sum(aa_dictionary[x]['C'] for x in clean) / L
    s = sum(aa_dictionary[x]['S'] for x in clean) / L
    mw = sum(aa_dictionary[x]['MW'] for x in clean) / L
    return round(n, 4), round(c, 4), round(s, 4), round(mw, 4)

# ==================== 3. 主循环 ====================

genome_dir = "${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_features/genomes"
protein_dir = "${PROJECT_ROOT:-/path/to/08ESKAPE}/46phage_features/proteins"

with open("Molecular_Final_Full_Stats.txt", "w") as out:
    out.write("Sample\tGC\tSize_bp\tNc\tGC_rank\tNARSC\tCARSC\tSARSC\tAvgMW\n")
    
    for f in os.listdir(genome_dir):
        if f.endswith((".fna", ".fa")):
            # 处理 DNA
            dna_recs = list(SeqIO.parse(os.path.join(genome_dir, f), "fasta"))
            dna_full = "".join(str(r.seq) for r in dna_recs)
            gc = round(gc_fraction(dna_full)*100, 4)
            size = len(dna_full)
            nc, rank = calculate_NC_and_Rank(dna_full)
            
            # 处理蛋白质
            p_f = f.replace(".fna", ".faa").replace(".fa", ".faa")
            p_path = os.path.join(protein_dir, p_f)
            if os.path.exists(p_path):
                aa_recs = list(SeqIO.parse(p_path, "fasta"))
                aa_full = "".join(str(r.seq) for r in aa_recs)
                narsc, carsc, sarsc, mw = calculate_protein_metrics(aa_full)
            else:
                narsc = carsc = sarsc = mw = "Nan"
                
            out.write(f"{f}\t{gc}\t{size}\t{nc}\t{rank}\t{narsc}\t{carsc}\t{sarsc}\t{mw}\n")
            print(f"Processed: {f}")

print("\n[Complete] 结果已保存至 Molecular_Final_Full_Stats.txt")
