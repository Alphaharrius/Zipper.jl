#!/bin/bash
#SBATCH --nodelist=node05
#SBATCH --job-name=GMERA/chern/systemsize128/local2_nn=0.1i  # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=64        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=10G         # memory per cpu-core (4G is default)
#SBATCH --time=72:00:00          # total run time limit (HH:MM:SS)

module purge
module load julia

julia --threads=auto localsize2.jl