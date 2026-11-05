
# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 01_genome_catalogue / 01_dereplication
# dRep dereplication of the 11 947 ESKAPE genomes
# Manuscript      : Fig. 1A
# Original path   : pkusz_CPU:2drep/0readme.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep/renamed_genome

cp ${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep/res_ESKAPE/dereplicated_genomes/*.fa ./
#!/bin/bash
for file in *.fa; do
    prefix=$(basename "$file" .fa)
    awk -v prefix="$prefix" '
    BEGIN {counter=1}
    /^>/ {
        if($0 ~ /plasmid/) {
            printf ">%s-%d-plasmid\n", prefix, counter++;
        } else {
            printf ">%s-%d\n", prefix, counter++;
        }
        next
    }
    {print}' "$file" > "${file%.fa}_renamed.fa"
done

cd ..

cat ${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep/renamed_genome/*renamed.fa > all_4538_ESKAPE_genomes.fa

sed '/^>/ s/\./_/g' all_4538_ESKAPE_genomes.fa > all_4538_ESKAPE_genomes_renamed.fa


cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep/renamed_4538genomes
for file in *.fa; do
sed -i '/^>/ s/\./_/g' ${file}

done
