#!/bin/bash
#SBATCH --account=PAS2493
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0 
#SBATCH -t 0-00:30:00
#SBATCH -J defModel
#SBATCH -o out.txt
#SBATCH -e err.txt

set -xe

module load oneapi/2024.0.2

export MV2_ENABLE_AFFINITY=0
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

echo "$SLURM_JOB_NUM_NODES nodes with $SLURM_CPUS_ON_NODE processors each." > log.txt

cd $SLURM_SUBMIT_DIR
date >> log.txt
make clean &>> log.txt
make -j$SLURM_CPUS_ON_NODE &>> log.txt

date >> log.txt
./bin/defaultModel &>> log.txt
date >> log.txt
echo "End." &>> ./log.txt

