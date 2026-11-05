
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 01_ice_integron_isescan
# Insertion sequence detection with ISEScan
# Manuscript      : Fig. 2A
# Original path   : pkusz_CPU:6ISEscan/00readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

find ${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep/renamed_4538genomes -type f -exec realpath {} \; > 4538MAG_path.txt
