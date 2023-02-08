#!/bin/bash
#SBATCH --job-name=RNASeq_nf
#SBATCH --cpus-per-task=100
# #SBATCH --mem-per-cpu=8G
#SBATCH --nodes=1
# #SBATCH --exclusive
#SBATCH --time=5-02:00:00
#SBATCH --output=R-_%j.log
#SBATCH --error=R-_%j.err

module load singularity-3.8.3-gcc-11.2.0-rlxj6fi



cd $SLURM_SUBMIT_DIR

#singularity exec /mnt/beegfs/idevillasante/apps/rocker/images/meth-dev.sif R CMD BATCH --no-save --no-restore '--args ncores=$SLURM_JOB_CPUS_PER_NODE  ' $dir/run.R slurm-batch-job.out


singularity exec /mnt/beegfs/idevillasante/apps/rocker/images/meth-dev.sif ./run.R -s multicore -n $SLURM_JOB_CPUS_PER_NODE  >& slurm-batch-job.out