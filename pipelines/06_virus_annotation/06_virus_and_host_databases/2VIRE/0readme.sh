# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 06_virus_and_host_databases
# VIRE download and ORF prediction for the ESKAPE virus references
# Manuscript      : Fig. 4B
# Original path   : lfat1t01:2virus_database/2VIRE/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE

grep 'Enterococcus faecium' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE/vire_predicted_host.tsv | wc -l
#0
grep 'Staphylococcus aureus' ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE/vire_predicted_host.tsv | wc -l
#464

#宿主预测信息与病毒contig ID无法匹配！

# 进入目录
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE/vire_genomes

# 使用find + xargs + gzip管道流式合并
find . -name "*.fna.gz" -type f -print0 | \
    xargs -0 -I {} sh -c 'gzip -dc "{}" && printf "\n"' | \
    cat > ../vire_all_genomes_combined.fasta
	
	
