# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 09_phylogeny / 02_special_genes
# Phylogenies of the special genes (e.g. VanX) found in the new viruses
# Manuscript      : Fig. 6
# Original path   : lfat1t01:61special_gene_phylogeny/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/61special_gene_phylogeny

#Dps蛋白是细菌用来应对氧化损伤和环境压力

#VanX 是一种锌依赖的 D-Ala-D-Ala 二肽酶，专门水解细胞壁合成前体中的 D-丙氨酰-D-丙氨酸（D-Ala-D-Ala）二肽

#Dps
seqkit grep -p 'vire_genome_12877_39' ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/ESKAPE_complete_phage.faa > phage_Dps.faa

#VanX
seqkit grep -p 'vire_genome_12877_30' ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/ESKAPE_complete_phage.faa > phage_VanX.faa


#Dps blast结果输出处理
query_len=149

awk -v L="$query_len" 'BEGIN{OFS="\t"; found=0} 
{
    # 提取字段（假设BLAST输出包含14列：标准12列 + staxids + stitle）
    qseqid = $1
    sseqid = $2
    pident = $3
    length_align = $4
    mismatch = $5
    gapopen = $6
    qstart = $7
    qend = $8
    sstart = $9
    send = $10
    evalue = $11
    bitscore = $12
    staxids = $13
    stitle = $14

    # 计算覆盖度（百分比）
    coverage = (qend - qstart + 1) / L * 100

    # 筛选条件
    if (pident > 30 && coverage > 60) {
        # 输出原行并附加覆盖度列
        print $0, coverage
        found = 1
    }
}
END {
    if (found == 0) {
        print "未找到满足 identity > 30% 且 coverage > 60% 的条目" > "/dev/stderr"
    }
}' "${PROJECT_ROOT:-/path/to/08ESKAPE}/61special_gene_phylogeny/Dps_vs_nr_blast.raw" > ./Dps_vs_nr_blast_2.tsv

cut -f2 ./Dps_vs_nr_blast_2.tsv | seqkit grep -f - ${DB_ROOT:-/path/to/databases}/NR/nr.faa > NR_extracted_Dps.faa

cat NR_extracted_Dps.faa phage_Dps.faa > merged_Dps.faa

#############################################################################
#VanX blast结果输出处理
seqkit stats ${PROJECT_ROOT:-/path/to/08ESKAPE}/61special_gene_phylogeny/phage_VanX.faa
query_len=268

awk -v L="$query_len" 'BEGIN{OFS="\t"; found=0} 
{
    # 提取字段（假设BLAST输出包含14列：标准12列 + staxids + stitle）
    qseqid = $1
    sseqid = $2
    pident = $3
    length_align = $4
    mismatch = $5
    gapopen = $6
    qstart = $7
    qend = $8
    sstart = $9
    send = $10
    evalue = $11
    bitscore = $12
    staxids = $13
    stitle = $14

    # 计算覆盖度（百分比）
    coverage = (qend - qstart + 1) / L * 100

    # 筛选条件
    if (pident > 30 && coverage > 60) {
        # 输出原行并附加覆盖度列
        print $0, coverage
        found = 1
    }
}
END {
    if (found == 0) {
        print "未找到满足 identity > 30% 且 coverage > 60% 的条目" > "/dev/stderr"
    }
}' "${PROJECT_ROOT:-/path/to/08ESKAPE}/61special_gene_phylogeny/VanX_vs_nr_blast2.raw" > ./VanX_vs_nr_blast_2.tsv

cut -f2 ./VanX_vs_nr_blast_2.tsv | seqkit grep -f - ${DB_ROOT:-/path/to/databases}/NR/nr.faa > NR_extracted_VanX.faa

cat NR_extracted_VanX.faa phage_VanX.faa > ./merged_VanX.faa

cat ./merged_VanX.faa rcsb_pdb_1R44.fasta rcsb_pdb_5HNM.fasta > ./merged_VanX2.faa

seqkit grep -p 'IMGVR_UViG_3300035126_000413|3300035126|Ga0335057_0000171_30' ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/ESKAPE_complete_phage.faa > ab_v_2_peptidase.faa

cat ./merged_VanX.faa ab_v_2_peptidase.faa Listeria_phageA118_peptidase.faa Listeria_phageA500_peptidase.faa > ./merged_VanX3.faa
