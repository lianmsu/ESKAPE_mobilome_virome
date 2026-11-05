
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 00_notes
# Integron detection with IntegronFinder
# Manuscript      : Fig. 2A
# Original path   : pkusz_CPU:4integronfinder/integron_filter/readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/4integronfinder/integron_filter

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/4integronfinder/res/Results_Integron_Finder_all_4538_ESKAPE_genomes_renamed/*.integrons ./

grep 'complete' all_4538_ESKAPE_genomes_renamed.integrons | cut -f2 | sort | uniq | wc -l
#1594

grep 'complete' all_4538_ESKAPE_genomes_renamed.integrons | cut -f2 | sort | uniq >> complete_integron_contig.id

head -n1 all_4538_ESKAPE_genomes_renamed.integrons > all_4538_ESKAPE_genomes_renamed_complete_integrons.tsv
grep 'complete' all_4538_ESKAPE_genomes_renamed.integrons >> all_4538_ESKAPE_genomes_renamed_complete_integrons.tsv
