# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 01_mge_crispr_spacers
# Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets)
# Manuscript      : Fig. 2B
# Original path   : lfat1t01:35MGE_spacer_competition/2blastn_db/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/2blastn_db

mamba activate blast2.16


#构建成功
makeblastdb -in ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/all_combined_ICEs.fa -dbtype nucl -out ./ICE.db -max_file_sz '4GB'

makeblastdb -in ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/all_plasmid_contig.fa -dbtype nucl -out ./Plasmid.db -max_file_sz '4GB'

makeblastdb -in ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/ESKAPE_complete_integrons_clean.fa -dbtype nucl -out ./Integron.db -max_file_sz '4GB'

makeblastdb -in ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/ESKAPE_genome_prophage_filtered.fa -dbtype nucl -out ./Prophage.db -max_file_sz '4GB'
