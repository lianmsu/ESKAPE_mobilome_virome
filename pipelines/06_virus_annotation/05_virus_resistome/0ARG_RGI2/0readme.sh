# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 05_virus_resistome
# ARG annotation of the viruses with RGI/CARD
# Manuscript      : Fig. 3F
# Original path   : lfat1t01:48virus_RGI/0ARG_RGI2/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/48virus_RGI/0ARG_RGI2
cut -f1,9,15,16,17,28 rgiargs.txt > rgiargs_detail_info.txt
