#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH -n 1
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -J all_phage_TDI
#SBATCH -o all_phage_TDI.out
#SBATCH -e all_phage_TDI.err

python compareTDI.py sequenced_phage_map_B3.txt ../../figures/TDI_individual_clusters/TDI_B3 "TDI first 10 phage of B3"
python compareTDI.py sequenced_phage_map_B4.txt ../../figures/TDI_individual_clusters/TDI_B4 "TDI first 10 phage of B4"
python compareTDI.py sequenced_phage_map_B5.txt ../../figures/TDI_individual_clusters/TDI_B5 "TDI first 10 phage of B5"
python compareTDI.py sequenced_phage_map_C1.txt ../../figures/TDI_individual_clusters/TDI_C1 "TDI first 10 phage of C1"
python compareTDI.py sequenced_phage_map_C2.txt ../../figures/TDI_individual_clusters/TDI_C2 "TDI first 10 phage of C2"
python compareTDI.py sequenced_phage_map_D1.txt ../../figures/TDI_individual_clusters/TDI_D1 "TDI first 10 phage of D1"
python compareTDI.py sequenced_phage_map_D2.txt ../../figures/TDI_individual_clusters/TDI_D2 "TDI first 10 phage of D2"
python compareTDI.py sequenced_phage_map_E.txt ../../figures/TDI_individual_clusters/TDI_E "TDI first 10 phage of E"
