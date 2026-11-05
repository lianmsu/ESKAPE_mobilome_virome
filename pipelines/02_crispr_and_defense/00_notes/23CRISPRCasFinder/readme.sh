# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 00_notes
# CRISPR-Cas detection with CRISPRCasFinder (level 4 arrays)
# Manuscript      : Fig. 2B
# Original path   : lfat1t01:23CRISPRCasFinder/readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/23CRISPRCasFinder

cut -f1 ${PROJECT_ROOT:-/path/to/08ESKAPE}/22minced/minced_res/all_4538_ESKAPE_genomes_minced_clean.gff \
    | sort \
    | uniq \
    | seqkit grep -f - \
    ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/all_4538_ESKAPE_genomes_renamed.fa \
    > minced_crispr_contig.fa


cp -r ${HOME_ROOT:-/path/to/home}/software/CRISPRCasFinder/CRISPRCasFinder-master ./
