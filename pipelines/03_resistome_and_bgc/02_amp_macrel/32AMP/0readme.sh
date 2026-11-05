
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 03_resistome_and_bgc / 02_amp_macrel
# Antimicrobial peptide prediction with Macrel
# Manuscript      : Fig. 1D
# Original path   : lfat1t01:32AMP/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/32AMP

#prodigal -p single -i ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/30minced_5455MAG/merged_5455MAG.fa -o 5455MAG.gff -a 5455MAGproteins.faa -d 5455MAGgenes.fna -g 11 -n -q
# 筛选长度 ≥15bp 且完整的 ORFs
awk '{if ($5-$4 >= 33 && $5-$4 <= 303 && $0 ~ /partial=00/) print}' ${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1/all_4538_ESKAPE_genomes_renamed.prt.gff > 4538ESKAPE_filtered_ORFs.gff

-p meta	启用宏基因组模式（适应复杂样本）
-g 11	使用宽松的起始密码子表（增加 sORF 检出）
-n	关闭动态长度阈值（避免自动过滤短 ORFs）
-c	仅保留完整 ORFs（起始+终止密码子）
-m	允许无起始密码子的 ORFs（需后续过滤）
-q	静默模式（不输出进度信息）

#根据文件提取文件5455MAGproteins.faa中对应的蛋白序列。
# 例如FASTA头格式：>CJ_14-10BD_SAB-bin.1-1_3
awk -F"\t" '{
    split($9, a, ";");
    for (i in a) {
        if (a[i] ~ /^ID=/) {
            sub(/^ID=/, "", a[i]);
			split(a[i], parts, "_");
            print $1 "_" parts[4];
            break;
        }
    }
}' 4538ESKAPE_filtered_ORFs.gff | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1/all_4538_ESKAPE_genomes_renamed.prt -o extracted_proteins.faa
