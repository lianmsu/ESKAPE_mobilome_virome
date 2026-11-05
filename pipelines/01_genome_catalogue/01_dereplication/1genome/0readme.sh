# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 01_genome_catalogue / 01_dereplication
# Download and rename the complete ESKAPE genomes from NCBI RefSeq
# Manuscript      : Fig. 1A
# Original path   : pkusz_CPU:1genome/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome

mkdir -p all_genomes

find ./ -type f -name "*.fna" | wc -l

find ./ -type f -name "*.fna" -exec sh -c '
  for f do
    base=$(basename "$f")          # 去除前缀
	echo $base
	echo $f
    cp "$f" "./all_genomes/$base"
  done
' sh {} +


