# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 10_experiments / 03_host_defense
# Defense system content of the experimental host strains
# Original path   : lfat1t01:70experiments/03host_DF/0_readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/03host_DF

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/70experiments/02host_CRBC_phylogeny/GCA* > ./merged_10_genome.fna

source activate defensefinder
defense-finder run -a ./merged_10_genome.fna --preserve-raw -o ./DF_res


