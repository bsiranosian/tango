#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH -n 1
#SBATCH --qos=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -J all_phage_TDI2
#SBATCH -o all_phage_TDI2.out
#SBATCH -e all_phage_TDI2.err

python compareTDI.py sequenced_phage_map_F1.txt ../../figures/TDI_individual_clusters/TDI_F1 "TDI first 10 phage of F1"
python compareTDI.py sequenced_phage_map_F2.txt ../../figures/TDI_individual_clusters/TDI_F2 "TDI first 10 phage of F2"
python compareTDI.py sequenced_phage_map_F3.txt ../../figures/TDI_individual_clusters/TDI_F3 "TDI first 10 phage of F3"
python compareTDI.py sequenced_phage_map_G.txt ../../figures/TDI_individual_clusters/TDI_G "TDI first 10 phage of G"
python compareTDI.py sequenced_phage_map_H1.txt ../../figures/TDI_individual_clusters/TDI_H1 "TDI first 10 phage of H1"
python compareTDI.py sequenced_phage_map_H2.txt ../../figures/TDI_individual_clusters/TDI_H2 "TDI first 10 phage of H2"
python compareTDI.py sequenced_phage_map_I1.txt ../../figures/TDI_individual_clusters/TDI_I1 "TDI first 10 phage of I1"
python compareTDI.py sequenced_phage_map_I2.txt ../../figures/TDI_individual_clusters/TDI_I2 "TDI first 10 phage of I2"
python compareTDI.py sequenced_phage_map_J.txt ../../figures/TDI_individual_clusters/TDI_J "TDI first 10 phage of J"
python compareTDI.py sequenced_phage_map_S.txt ../../figures/TDI_individual_clusters/TDI_S "TDI first 10 phage of S"
python compareTDI.py sequenced_phage_map_Singleton.txt ../../figures/TDI_individual_clusters/TDI_Singleton "TDI first 10 phage of Singleton"
python compareTDI.py sequenced_phage_map_T.txt ../../figures/TDI_individual_clusters/TDI_T "TDI first 10 phage of T"