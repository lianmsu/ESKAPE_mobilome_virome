#!/bin/bash
#SBATCH -o job.gtdbtk.%j.%N.out 
#SBATCH --partition=fat
#SBATCH -J GTDB_crbc
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=500GB

#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 10_experiments / 05_extra_validation
# Additional GTDB-Tk validation of the experimental strains
# Original path   : pkusz_CPU:70_experiments/dl_GTDB_classify.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate gtdbtk-2.4.0

binning_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/70_experiments/CRBC
#out_path=${HOME_ROOT:-/path/to/home}/11basic_lake/02DB/07GTDB



start_time=$(date +%s)
cd ${PROJECT_ROOT:-/path/to/08ESKAPE}/70_experiments
mkdir -p res_EcoliMG1655
gtdbtk classify_wf --cpus $SLURM_NTASKS --pplacer_cpus  $SLURM_NTASKS  -x fna \
    --genome_dir ${binning_path}/  \
	--out_dir ./res_EcoliMG1655 --skip_ani_screen
    #--mash_db ${HOME_ROOT:-/path/to/home}/02database/gtdbtk_release214/mash  \
    

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "running time is $(($cost_time/60)) min" >> run_time_gtdb_${SLURM_JOBID}.log



#mkdir ${out_path}/res

#gtdbtk classify_wf --cpus 14 \
#--genome_dir ${binning_path}/res_merged/dereplicated_genomes \
#--out_dir ${out_path}/res \
#-x fa \
#--pplacer_cpus 1

