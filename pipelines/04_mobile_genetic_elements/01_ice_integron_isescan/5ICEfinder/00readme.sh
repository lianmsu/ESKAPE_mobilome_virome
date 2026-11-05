
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 01_ice_integron_isescan
# ICE detection with ICEfinder2 (batches of 1000 and the full set)
# Manuscript      : Fig. 2A
# Original path   : pkusz_CPU:5ICEfinder/00readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/5ICEfinder
awk -F "," 'NR>1{print $1}' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep/res_ESKAPE/data_tables/Widb.csv > 4538MAG.id
