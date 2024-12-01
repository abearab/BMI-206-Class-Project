#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -j y               # STDERR and STDOUT should be joined
#$ -l mem_free=2G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=5G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=4:00:00   # job requires up to 24 hours of runtime
#$ -t 1-1000          # array job with 10 tasks (remove first '#' to enable)
#$ -r y               # if job crashes, it should be restarted

# Parse command-line arguments
export num_cores=6         # Number of cores
export file_SCENT_obj="./SCENT_obj_all.rds"    # Path to SCENT object
export celltype="Tcell"    # Cell type
export bin="TRUE"
export regr="poisson"
export output_dir="./"     # Output directory

module load CBI
module load r
Rscript SCENT_parallelization.R $SGE_TASK_ID ${num_cores} ${file_SCENT_obj} ${celltype} ${regr} ${bin} ${output_dir}
