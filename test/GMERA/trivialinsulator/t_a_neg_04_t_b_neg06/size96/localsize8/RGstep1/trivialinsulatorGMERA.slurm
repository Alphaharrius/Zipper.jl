#!/bin/bash
#SBATCH --job-name=GMERA/trivial/t_a_neg_04_t_b_neg_06/size96/localsize8/RGstep1    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=8        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=50G         # memory per cpu-core (4G is default)
#SBATCH --time=72:00:00          # total run time limit (HH:MM:SS)

module purge
module load julia

julia --threads=auto trivialinsulatorGMERA.jl 
