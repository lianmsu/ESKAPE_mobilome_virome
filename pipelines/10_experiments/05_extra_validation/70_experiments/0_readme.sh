# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 10_experiments / 05_extra_validation
# Additional GTDB-Tk validation of the experimental strains
# Original path   : pkusz_CPU:70_experiments/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/70_experiments

#从本部服务器复制过来10个ESKAPE基因组，与CRBC基因组构建系统发育树

#对CRBC基因组进行软链接
ln -s ${HOME_ROOT:-/path/to/home}/Culturome/03_CRBC/01data/CRBC_6699genome/* ${PROJECT_ROOT:-/path/to/08ESKAPE}/70_experiments/CRBC

mv ./GCA* ./CRBC

###################
#手动放入一个大肠杆菌的基因组，进行建树
