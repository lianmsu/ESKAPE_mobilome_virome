# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 05_virus_discovery / 03_clustering
# CheckV + ANI clustering of the combined virus set
# Manuscript      : Fig. 3A
# Original path   : lfat1t01:44virus_checkv_cluster/0checkv_cluster.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/44virus_checkv_cluster

makeblastdb -in ${PROJECT_ROOT:-/path/to/08ESKAPE}/43combined_all_virus/all_ESKAPE_phage.fa -dbtype nucl -out ./my_db
blastn -query ${PROJECT_ROOT:-/path/to/08ESKAPE}/43combined_all_virus/all_ESKAPE_phage.fa -db ./my_db -outfmt '6 std qlen slen' -max_target_seqs 10000 -out my_blast.tsv -num_threads 10
./anicalc.py -i ./my_blast.tsv -o ./my_ani.tsv
./aniclust.py --fna ${PROJECT_ROOT:-/path/to/08ESKAPE}/43combined_all_virus/all_ESKAPE_phage.fa --ani ./my_ani.tsv --out ./my_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0

cut -f1 ${PROJECT_ROOT:-/path/to/08ESKAPE}/44virus_checkv_cluster/my_clusters.tsv | seqkit grep -f - ${PROJECT_ROOT:-/path/to/08ESKAPE}/43combined_all_virus/all_ESKAPE_phage.fa > ./checkv_clustered_phage.fa

sed '/^>/ s/provirus/prophage/g' checkv_clustered_phage.fa > checkv_clustered_modified_phage.fa
