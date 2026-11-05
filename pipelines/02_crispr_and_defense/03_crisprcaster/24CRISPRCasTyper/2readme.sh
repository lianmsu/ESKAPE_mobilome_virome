# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 02_crispr_and_defense / 03_crisprcaster
# CRISPR-Cas subtyping with CRISPRCasTyper
# Manuscript      : Fig. 1C
# Original path   : lfat1t01:24CRISPRCasTyper/2readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/24CRISPRCasTyper
cut -f2 selected_CRISPRCas.txt | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/24CRISPRCasTyper/all_CRISPR_level4.fa > selected_plot_contig.fa


source activate cctyper

cctyper selected_plot_contig.fa ./res_selected_plot -t 10 --prodigal single --no_grid

