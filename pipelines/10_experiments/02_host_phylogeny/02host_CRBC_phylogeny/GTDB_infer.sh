#!/bin/bash
#SBATCH -o job.gtdbtk-infer.%j.%N.out 
#SBATCH --partition=cpu
#SBATCH -J infer
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28

#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 10_experiments / 02_host_phylogeny
# Phylogeny of the CRBC host strains (GTDB-Tk)
# Original path   : lfat1t01:70experiments/02host_CRBC_phylogeny/GTDB_infer.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate gtdbtk-2.5.2

cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/1genome/4GTDB

start_time=$(date +%s)

gtdbtk infer --msa_file gtdbtk.bac120.user_msa.fasta.gz --out_dir ./infer --prefix 4538ESKAPE --cpus 28

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "running time is $(($cost_time/60)) min" >> run_time_gtdb_${SLURM_JOBID}.log

