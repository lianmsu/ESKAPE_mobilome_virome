
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 00_notes
# Defense system annotation of the ESKAPE genomes with DefenseFinder
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:30Defensefinder/res1/readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1

sed -i 's/\*//g' all_4538_ESKAPE_genomes_renamed.prt

grep '>' all_4538_ESKAPE_genomes_renamed.prt | \
awk -F ' # ' 'BEGIN {OFS="\t"; print "##gff-version 3"} {
    # 提取seqid（去掉">"和末尾的"_数字"）
    seqid = substr($1, 2);
    sub(/_[0-9]+$/, "", seqid);
    
    # 提取protein_id（最后的"_数字"）
    protein_id = $1;
    sub(/^>/, "", protein_id);
    
    # 处理链方向
    strand = ($4 == "1" ? "+" : "-");
    
    # 构造attributes
    attrs = "ID=" protein_id ";";
    n = split($5, a, ";");
    for (i=1; i<=n; i++) {
        if (a[i] ~ /^ID=/) continue;  # 跳过重复的ID
        attrs = attrs a[i] ";";
    }
    
    # 输出GFF行
    print seqid, "DefenseFinder", "CDS", $2, $3, ".", strand, ".", attrs;
}' > all_4538_ESKAPE_genomes_renamed.prt.gff

#根据文件${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/res1/merged_5467MAG_defense_finder_genes.tsv的hit_id列从文件${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/res1/merged_5467MAG.prt_defensefinder.prt中提取对应的蛋白序列到新的文件：All_DG.faa
awk -F"\t" 'NR>1 {print $2}' ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/res1/merged_5467MAG_defense_finder_genes.tsv > hit_ids.txt
seqkit grep -f hit_ids.txt ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs/res1/merged_5467MAG.prt_defensefinder.prt > All_DG.faa
seqkit stat All_DG.faa
