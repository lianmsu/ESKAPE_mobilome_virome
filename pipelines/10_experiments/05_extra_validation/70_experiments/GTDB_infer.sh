#!/bin/bash
#SBATCH -o job.gtdbtk.%j.%N.out 
#SBATCH --partition=cpu
#SBATCH -J GTDB_infer_27578
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=500GB

#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 10_experiments / 05_extra_validation
# Additional GTDB-Tk validation of the experimental strains
# Original path   : pkusz_CPU:70_experiments/GTDB_infer.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate gtdbtk-2.4.0



start_time=$(date +%s)
cd .
gtdbtk infer --msa_file gtdbtk.bac120.user_msa.fasta.gz --out_dir ./infer --prefix bin27578 --cpus 48

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "running time is $(($cost_time/60)) min" >> run_time_gtdb_${SLURM_JOBID}.log



#mkdir ${out_path}/res

#gtdbtk classify_wf --cpus 14 \
#--genome_dir ${binning_path}/res_merged/dereplicated_genomes \
#--out_dir ${out_path}/res \
#-x fa \
#--pplacer_cpus 1

