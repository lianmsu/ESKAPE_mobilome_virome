#!/bin/sh
#SBATCH -o job.tree.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J vanx3
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28


# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 09_phylogeny / 02_special_genes
# Phylogenies of the special genes (e.g. VanX) found in the new viruses
# Manuscript      : Fig. 6
# Original path   : lfat1t01:61special_gene_phylogeny/2_iqtree_VanX3.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/61special_gene_phylogeny
source activate iqtree2


linsi --thread 25 ./merged_VanX3.faa > merged_VanX3_mafft.afa
clipkit merged_VanX3_mafft.afa -m gappy -g 0.9
iqtree2 -s merged_VanX3_mafft.afa.clipkit -m MFP+MERGE -B 1000 -nt AUTO --nmax 4000
