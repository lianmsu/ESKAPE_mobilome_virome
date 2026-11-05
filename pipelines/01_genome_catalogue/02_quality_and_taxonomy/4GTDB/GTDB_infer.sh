#!/bin/bash
#SBATCH -o job.gtdbtk-infer.%j.%N.out 
#SBATCH --partition=cpu
#SBATCH -J infer
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28

#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 01_genome_catalogue / 02_quality_and_taxonomy
# GTDB-Tk taxonomy and reference tree of the dereplicated genomes
# Manuscript      : Fig. 1B
# Original path   : lfat1t01:1genome/4GTDB/GTDB_infer.sh
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

