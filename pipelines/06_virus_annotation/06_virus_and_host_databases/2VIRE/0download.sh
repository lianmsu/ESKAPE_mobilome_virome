
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 06_virus_annotation / 06_virus_and_host_databases
# VIRE download and ORF prediction for the ESKAPE virus references
# Manuscript      : Fig. 4B
# Original path   : lfat1t01:2virus_database/2VIRE/0download.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

wget --no-check-certificate -O vire_metadata_env_all_samples.tsv 'https://spire.embl.de/vire/api/download_metadata/?environment=all_samples&data_type=metadata&format=txt'

wget --no-check-certificate -O vire_all_genomes.list 'https://spire.embl.de/vire/api/download_files/?environment=all_samples&file_type=genomes&format=txt'

# 1. 使用 wget 的 -i 参数从文件读取URL
wget -i ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE/vire_all_genomes.list

# 2. 上述命令会将文件保存为链接中的原始文件名。
#    如果要下载到指定目录并显示进度，可以：
wget -i ${PROJECT_ROOT:-/path/to/08ESKAPE}/2virus_database/2VIRE/vire_all_genomes.list -P ./vire_genomes/ --progress=bar
