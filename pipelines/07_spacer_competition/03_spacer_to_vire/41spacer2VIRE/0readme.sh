
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 03_spacer_to_vire
# BLAST of the host spacers against the VIRE database
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:41spacer2VIRE/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/41spacer2VIRE

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/25level4_minced/minced_res/level4_minced_spacers.fa ./




#spacer通过blast比对到phage 和plasmid

#筛选被spacer靶向的phage和plasmid（含有anti defense基因的）

cut -f1,2 ${PROJECT_ROOT:-/path/to/08ESKAPE}/41spacer2VIRE/spacers_VIRE.blastn.1mis95cov.out > spacers_VIRE.blastn.1mis95cov.out2
