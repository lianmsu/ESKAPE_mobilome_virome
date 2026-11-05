#!/bin/bash
#SBATCH -o job.gtdbtk.%j.%N.out 
#SBATCH --partition=fat
#SBATCH -J infer_4538
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=400GB

#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 01_genome_catalogue / 02_quality_and_taxonomy
# GTDB-Tk classification of the dereplicated genomes
# Manuscript      : Fig. 1B
# Original path   : pkusz_CPU:3gtdbtk/align/GTDB_infer.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate gtdbtk-2.4.0


cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/3gtdbtk/align
start_time=$(date +%s)
cd .
gtdbtk infer --msa_file gtdbtk.bac120.user_msa.fasta.gz --out_dir ./infer --prefix bin4538 --cpus 40

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "running time is $(($cost_time/60)) min" >> run_time_gtdb_${SLURM_JOBID}.log



#mkdir ${out_path}/res

#gtdbtk classify_wf --cpus 14 \
#--genome_dir ${binning_path}/res_merged/dereplicated_genomes \
#--out_dir ${out_path}/res \
#-x fa \
#--pplacer_cpus 1

