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
#SBATCH --job-name=hyp_vel_testruns_cpu
#SBATCH --time=16:00:00
#SBATCH --array=73-84
#SBATCH --gres=gpu:0
#SBATCH --gpus=0
echo "Program start: $(date)"

id=$(echo "$SLURM_ARRAY_TASK_ID")
source /home/lguelzow/hyp_vel/bin/activate

mkdir Results_$id
cd Results_$id
python3 -u "/home/lguelzow/half-mass/RNG-IC-List.py"
time python3 -u "/home/lguelzow/half-mass/HVS-SIM-paper-final.py"

echo "Program end: $(date)"
