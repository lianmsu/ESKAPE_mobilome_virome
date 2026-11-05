#!/bin/bash
#SBATCH -o job.gtdbtk.%j.%N.out 
#SBATCH --partition=cpu
#SBATCH -J gtdb
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=500GB

#SBATCH --mail-type=end

# ===========================================================================
# ESKAPE_mobilome_virome -- pipelines / 01_genome_catalogue / 02_quality_and_taxonomy
# GTDB-Tk classification of the dereplicated genomes
# Manuscript      : Fig. 1B
# Original path   : pkusz_CPU:3gtdbtk/01GTDB_classify.sh
#
# Site-specific paths were replaced by ${PROJECT_ROOT} (the 08ESKAPE
# directory), ${HOME_ROOT} (home of the companion projects), ${SHARED_HOME}
# and ${DB_ROOT} (reference databases); export them before submitting.
# ===========================================================================

source activate gtdbtk-2.4.0


#GEM
binning_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/2drep/renamed_4538genomes
out_path=${PROJECT_ROOT:-/path/to/08ESKAPE}/3gtdbtk

start_time=$(date +%s)
cd ${out_path}
#mkdir -p res
gtdbtk classify_wf --cpus $SLURM_NTASKS --pplacer_cpus $SLURM_NTASKS  -x fa \
    --genome_dir ${binning_path}/  \
	--out_dir ${out_path} --skip_ani_screen    

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "running time is $(($cost_time/60)) min" >> run_time_gtdb_${SLURM_JOBID}.log



