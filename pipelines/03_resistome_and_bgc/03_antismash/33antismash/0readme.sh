# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 03_resistome_and_bgc / 03_antismash
# Biosynthetic gene cluster prediction with antiSMASH
# Manuscript      : Fig. 1D
# Original path   : lfat1t01:33antismash/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

#把路径下所有以.gbk结尾的文件cat到一起。
#${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/11prokka_5467MAGs/res
find ${HOME_ROOT:-/path/to/home}/04Qinghai_Tibet_Plateau/11prokka_5467MAGs/res -name "*.gbk" -type f -exec cat {} + > 5467MAGcombined.gbk

