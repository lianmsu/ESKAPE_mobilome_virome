# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 02_spacer_to_imgvr
# BLAST of the host spacers against the IMG/VR database
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:40spacer2IMGVR/4res_filter.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/40spacer2IMGVR

# cp ${HOME_ROOT:-/path/to/home}/06mobilome_predators/62_NR_filtering/blastp_res_tax/contig_dominant_taxa_manuscript_filter.tsv ./

# #过滤文件spacers_virus.blastn.1mis95cov.out，只保留该文件中包含文件contig_dominant_taxa_manuscript_filter.tsv第一列字符串的行，不用严格匹配。
# # 提取contig_dominant_taxa_manuscript_filter.tsv的第一列内容
# cut -f1 contig_dominant_taxa_manuscript_filter.tsv | tail -n +2 > filter_patterns.txt

# # 使用grep进行包含匹配
# grep -f filter_patterns.txt spacers_virus.blastn.1mis95cov.out > filtered_spacers_virus.blastn.1mis95cov.out

# # 使用grep进行包含匹配
# grep -f filter_patterns.txt spacers_plasmid.blastn.1mis95cov.out > filtered_spacers_plasmid.blastn.1mis95cov.out

cut -f2 spacers_IMGVR.blastn.1mis95cov.out | sort | uniq | seqkit grep -f - ${HOME_ROOT:-/path/to/home}/06mobilome_predators/2IMG_VR_4/IMGVR_all_nucleotides.fna > spacers_IMGVR.blastn.1mis95cov.fa


cat ./spacer2IMGVR/spacers_IMGVR.blastn.out | awk '{if($5<=1 && $9>=100) print $0}' > spacers_IMGVR.blastn.1mis100cov.out

cut -f2 spacers_IMGVR.blastn.1mis100cov.out | sort | uniq | seqkit grep -f - ${HOME_ROOT:-/path/to/home}/06mobilome_predators/2IMG_VR_4/IMGVR_all_nucleotides.fna > spacers_IMGVR.blastn.1mis100cov.fa

cut -f1,2 ${PROJECT_ROOT:-/path/to/08ESKAPE}/40spacer2IMGVR/spacers_IMGVR.blastn.1mis100cov.out > spacers_IMGVR.blastn.1mis100cov.out2
