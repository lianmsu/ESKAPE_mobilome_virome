# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 02_special_mge
# ICEfinder2 on the selected plasmid/prophage set
# Manuscript      : Fig. 2E-H
# Original path   : pkusz_CPU:5_ICEfinder_special_MGE/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/5_ICEfinder_special_MGE

source activate icefinder

cp -r ${HOME_ROOT:-/path/to/home}/05plateau/21icefinder_5455MAGs/ICEfinder2_moban ./ICEfinder2_mge
ICE_workdir=./ICEfinder2_mge
cd ./ICEfinder2_mge

python ICEfinder2.py -i ${PROJECT_ROOT:-/path/to/08ESKAPE}/5_ICEfinder_special_MGE/selected_plasmid_prophage.fa -t Metagenome

