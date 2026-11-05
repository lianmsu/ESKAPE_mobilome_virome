# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 04_mobile_genetic_elements / 02_special_mge
# Selection of the special MGEs used in the case studies
# Manuscript      : Fig. 2E-H
# Original path   : lfat1t01:34special_MGE/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/34special_MGE
cat manual_selected_MGE.id | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/21genomad/res/all_4538_ESKAPE_genomes_renamed_summary/all_4538_ESKAPE_genomes_renamed_plasmid.fna > selected_plasmid.fa

cat manual_selected_MGE.id | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/21genomad/res/all_4538_ESKAPE_genomes_renamed_summary/all_4538_ESKAPE_genomes_renamed_virus.fna > selected_prophage.fa

cat selected_plasmid.fa selected_prophage.fa > selected_plasmid_prophage.fa

mamba activate prokka
prodigal -p meta -i selected_plasmid_prophage.fa \
		-f gff -o ./selected_plasmid_prophage.gff \
		-m -d ./selected_plasmid_prophage.fna \
		-a ./selected_plasmid_prophage.faa


source activate genomad
genomad_db=${DB_ROOT:-/path/to/databases}/genomad_db
#########################################################################
contig=${PROJECT_ROOT:-/path/to/08ESKAPE}/34special_MGE/selected_plasmid_prophage.fa
workdir=${PROJECT_ROOT:-/path/to/08ESKAPE}/34special_MGE

genomad end-to-end --cleanup ${contig} ${workdir}/genomad_res ${genomad_db}
