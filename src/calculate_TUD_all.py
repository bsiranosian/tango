#calculate TUD for all 663 phages as of 2014-04-20
from tetranucleotideAnalysis import *

#phageTUDCalc('../data/phagesDB/sequenced_phage_map.txt', '../data/all_phages_TUD_3.tsv', 3)
phageTUDCalc('../data/phagesDB/sequenced_phage_map.txt', '../data/all_phages_TUD_4.tsv', 4, RC=True)
#phageTUDCalc('../data/phagesDB/sequenced_phage_map.txt', '../data/all_phages_TUD_5.tsv', 5)
