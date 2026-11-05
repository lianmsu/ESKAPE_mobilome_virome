# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 01_minced
# CRISPR spacer extraction with minced on the ESKAPE genomes
# Manuscript      : Fig. 2B
# Original path   : lfat1t01:22minced/minced_res/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/22minced/minced_res
grep 'minced' ./all_4538_ESKAPE_genomes_minced.gff > ./all_4538_ESKAPE_genomes_minced_clean.gff
