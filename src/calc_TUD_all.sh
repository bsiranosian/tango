#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH -n 1
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -J all_phage_TUD
#SBATCH -o all_phage_TUD.out
#SBATCH -e all_phage_TUD.err

python calculate_TUD_all.py
