
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 07_spacer_competition / 02_spacer_to_imgvr
# BLAST of the host spacers against the IMG/VR database
# Manuscript      : Fig. 4C
# Original path   : lfat1t01:40spacer2IMGVR/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/40spacer2IMGVR

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/25level4_minced/minced_res/level4_minced_spacers.fa ./




#spacer通过blast比对到phage 和plasmid

#筛选被spacer靶向的phage和plasmid（含有anti defense基因的）

cut -f1,2 ${PROJECT_ROOT:-/path/to/08ESKAPE}/40spacer2IMGVR/spacers_IMGVR.blastn.1mis95cov.out > spacers_IMGVR.blastn.1mis95cov.out2
