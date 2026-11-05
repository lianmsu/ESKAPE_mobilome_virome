
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 03_resistome_and_bgc / 01_resistome_annotation
# ARG annotation: RGI/CARD and diamond against SARG
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:31ARG/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/31ARG
ln -s ${PROJECT_ROOT:-/path/to/08ESKAPE}/30Defensefinder/res1/all_4538_ESKAPE_genomes_renamed.prt ./4538genomes_prodigal.faa

#sed -i 's/\*//g' 4538genomes_defensefinder.faa

