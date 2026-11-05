# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 01_mge_crispr_spacers
# Spacer-based competition between MGE types (CRISPRCasTyper on MGE sets)
# Manuscript      : Fig. 2B
# Original path   : lfat1t01:35MGE_spacer_competition/1MGE_spacer/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/1MGE_spacer

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/CCTyper_res_prophage/spacers/*.fa > ./Prophage_spacers.fa

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/CCTyper_res_plasmid/spacers/*.fa > ./Plasmid_spacers.fa

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/CCTyper_res_integron/spacers/*.fa > ./Integron_spacers.fa

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/35MGE_spacer_competition/CCTyper_res_ICE/spacers/*.fa > ./ICE_spacers.fa



