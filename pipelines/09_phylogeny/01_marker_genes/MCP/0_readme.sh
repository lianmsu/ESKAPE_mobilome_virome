# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 09_phylogeny / 01_marker_genes
# hhsuite2 annotation and IQ-TREE phylogenies of MCP / portal / TerL
# Manuscript      : Fig. 6
# Original path   : lfat1t01:60marker_gene_phylogeny/MCP/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/MCP

tail -n+2 ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/MCP/filtered_32phage_MCP.tsv | cut -f6 | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/32_phage.faa > ./32phage_MCP.faa

#其他7个phage的MCP手动整理来
tail -n+2 ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/MCP/selected_7_phage_MCP.tsv | cut -f6 | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/49selected_virus_annotation/selected_phage.faa > ./7phage_MCP.faa

cat ./32phage_MCP.faa ./7phage_MCP.faa > 39phage_MCP.faa

#########################################################
tail -n+2 ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/MCP/filtered_phage_MCP.tsv | cut -f6 | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/60marker_gene_phylogeny/ESKAPE_complete_phage.faa > ./ESKAPE_complete_phage_MCP.faa
