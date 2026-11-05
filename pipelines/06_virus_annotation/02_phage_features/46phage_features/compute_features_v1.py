# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 02_phage_features
# Molecular features of the phages: GC, NARSC, length, quality
# Manuscript      : Fig. 3
# Original path   : lfat1t01:46phage_features/compute_features_v1.py
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

protein_dir = "proteins/"
output_AAU = []
output_pI = []

bins = np.arange(1, 14.25, 0.25)
bin_centers = bins[:-1] + 0.125

for file in os.listdir(protein_dir):

    if not file.endswith(".faa"):
        continue

    genome_id = file.replace(".faa", "")
    aa_counts = dict.fromkeys(list("ACDEFGHIKLMNPQRSTVWY"), 0)
    pI_values = []

    for record in SeqIO.parse(os.path.join(protein_dir, file), "fasta"):
        seq = str(record.seq)

        # AAU 统计
        for aa in aa_counts:
            aa_counts[aa] += seq.count(aa)

        # pI 计算
        analysed = ProteinAnalysis(seq)
        pI_values.append(analysed.isoelectric_point())

    total_aa = sum(aa_counts.values())
    aa_freq = {aa: aa_counts[aa]/total_aa for aa in aa_counts}
    aa_freq["Genome"] = genome_id
    output_AAU.append(aa_freq)

    # pI distribution
    hist, _ = np.histogram(pI_values, bins=bins)
    hist = hist / sum(hist)

    pI_dict = {"Genome": genome_id}
    for i, center in enumerate(bin_centers):
        pI_dict[f"pI_{round(center,2)}"] = hist[i]

    output_pI.append(pI_dict)

# 保存结果
pd.DataFrame(output_AAU).to_csv("AAU_matrix.csv", index=False)
pd.DataFrame(output_pI).to_csv("pI_distribution_matrix.csv", index=False)
