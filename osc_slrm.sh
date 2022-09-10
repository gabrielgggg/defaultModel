#!/bin/bash
#SBATCH --account=PAS2282
#SBATCH -p debug-40core
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0 
#SBATCH -t 0-01:00:00
#SBATCH -J defModel
#SBATCH -o out.txt
#SBATCH -e err.txt

set -xe

module load intel/2021.3.0

export MV2_ENABLE_AFFINITY=0
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

echo "$SLURM_JOB_NUM_NODES nodes with $SLURM_CPUS_ON_NODE processors each." > log.txt

cd $SLURM_SUBMIT_DIR
date >> log.txt
make clean &>> log.txt
make -j$SLURM_CPUS_ON_NODE &>> log.txt

date >> log.txt
# srun -n $SLURM_JOB_NUM_NODES -c $SLURM_CPUS_ON_NODE ./bin/defaultModel &>> log.txt
./bin/defaultModel &>> log.txt

date >> log.txt
echo "End." &>> ./log.txt

