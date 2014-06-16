#parameters 
import argparse

parser = argparse.ArgumentParser(description="This script calculates Tetranucleotide Usage Deviation (TUD) for multiple phage genomes. Originally it could only do 4-mers, but has now been generalized to all nucleotide lengths. The defualt function of this script is to save a nexus distance file that can be used for tree buliding with splitstree and other programs. The only required inputs are the location of a configuration file and location to save the resulting nexus distance file. The configuration file is a comma separated file of phage info with 2 necessary fields and 2 optional fields: \n name,fastaPath,subsetStart,subsetEnd\nEach phage to be compared should be separated by a new line. Additional arguments provide more control over the configuration the underlying calculations. If you want to save a csv file of the resulting data, specify one with the --s option. A distance matrix can be saved by specifying a file name with the --d option.")
parser.add_argument("nameFile", action="store", metavar="nameFile", help="The file name of the comma separated file defining phage information")
parser.add_argument("saveName", action="store", metavar="saveName", help="The file to save the resulting nexus to.")
parser.add_argument("--subset", action="store", metavar="subset", default="False", help="set to True to plot only calculate for the regions defined in the input file. If True, each line must have integers in fields 3 and 4 that represent the genomic region of each phage to compare. This allows you to subset on regions containing genes, etc for each phage.")
parser.add_argument("--k", action="store", metavar="k", default="4", help="Can also use this to compute 2,3,5-mers, etc.")
parser.add_argument("--s", action="store", metavar="saveFile", default=None, help="specify a filename to save resulting usage deviation data to. Be careful because these files can get big for higher values of k!")
parser.add_argument("--d", action="store", metavar="distanceFile", default=None, help="specify a filename to save resulting distance matrix to")
args=parser.parse_args()

nameFile=args.nameFile
saveName=args.saveName
if args.subset == "True": subset=True
else: subset=False
k=int(args.k)
saveFile = args.saveFile
saveFile = args.distanceFile

from tetranucelotideAnalysis import 