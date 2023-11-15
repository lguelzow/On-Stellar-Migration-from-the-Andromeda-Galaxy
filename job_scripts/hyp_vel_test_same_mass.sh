#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --cores-per-socket=1
#SBATCH --tasks-per-node=1
#SBATCH --output=log/out%a
#SBATCH --error=log/err%a
#SBATCH --mail-type=all
#SBATCH --mail-user=lguelzow@physik.uni-bielefeld.de
#SBATCH --job-name=hyp_vel_same-mass
#SBATCH --time=16:00:00
#SBATCH --array=21-40
#SBATCH --gres=gpu:0
#SBATCH --gpus=0
echo "Program start: $(date)"

id=$(echo "$SLURM_ARRAY_TASK_ID")
source /home/lguelzow/hyp_vel/bin/activate

mkdir Results_$id
cd Results_$id
python3 -u "/home/lguelzow/paper_update/On-Stellar-Migration-from-the-Andromeda-Galaxy/simulation-code/RNG_initial_conditions.py"
time python3 -u "/home/lguelzow/paper_update/On-Stellar-Migration-from-the-Andromeda-Galaxy/simulation-code/HVS_trajectory_simulation.py"

echo "Program end: $(date)"
