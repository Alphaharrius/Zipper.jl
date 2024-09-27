#!/bin/bash
#SBATCH --job-name=GMERA/chern/systemsize96/nnhoping_0.4i  # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=64        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=22G         # memory per cpu-core (4G is default)
#SBATCH --time=336:00:00          # total run time limit (HH:MM:SS)

module purge
module load julia

julia --threads=auto localsizes_implement.jl 