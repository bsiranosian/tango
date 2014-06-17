#parameters 
import argparse

parser = argparse.ArgumentParser(description="This script produces a Tetranucleotide Difference Index (TDI) plot for multiple phage genomes. Originally it could only do 4-mers, but has now been generalized to all nucleotide lengths. The only required inputs are the location of a configuration file, location to save the resulting figure and a title for the plot. The configuration file is a comma separated file of phage info with 2 necessary fields and 2 optional fields: \n name,fastaPath,subsetStart,subsetEnd\nEach phage to be compared should be separated by a new line. Additional arguments provide more control over the configuration of the plot and the underlying calculations. If you want to save a csv file of the resulting data, specify one with the --saveFile option. The length of kmer used can be changed with the --k option.")
parser.add_argument("nameFile", action="store", metavar="nameFile", help="The file name of the comma separated file defining phage information")
parser.add_argument("saveName", action="store", metavar="saveName", help="The file to save the resulting image to.")
parser.add_argument("title", action="store", metavar="title", help="The title for the resulting plot.")
parser.add_argument("--maxNum", action="store", metavar="maxnum", default="10", help="maximum number of phage to plot on the same figure. first maxnum of input file will be chosen. default: 10")
parser.add_argument("--xScale", action="store", metavar="xScale", default="False", help="Set to True to scale the x axis of the plot to relative genome position.")
parser.add_argument("--subset", action="store", metavar="subset", default="False", help="set to True to plot only the regions defined in the input file. If True, each line must have integers in fields 3 and 4 that represent the genomic region of each phage to compare")
parser.add_argument("--windowSize", action="store", metavar="windowSize", default="5000", help="The size of the window to compute TDI within. Default of 5000bp is used in Pride et. al 2006")
parser.add_argument("--stepSize", action="store", metavar="stepSize", default="1000", help="How many bases to move the window along the genome at each iteration. Default of 1000bp is used in Pride et. al 2006")
parser.add_argument("--k", action="store", metavar="k", default="4", help="Can also use this to compute 2,3,5-mers, etc.")
parser.add_argument("--saveFile", action="store", metavar="saveFile", default=None, help="specify a filename to save resulting data to")
args=parser.parse_args()

nameFile=args.nameFile
saveName=args.saveName
title=args.title
maxNum=int(args.maxNum)
if args.xScale == "True": xScale=True
else: xScale=False
if args.subset == "True": subset=True
else: subset=False
windowSize=int(args.windowSize)
stepSize=int(args.stepSize)
k=int(args.k)
saveFile = args.saveFile


from tetranucleotideAnalysis import *
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import numpy as np
import csv
#compareTDI(): plots multiple tetranucleotide difference index signals on the same plot. 
# takes as input a  comma separated file of phage info with 2 necessary fields and 2 optional fields:
# 	name ,fastaPath,subsetStart,subsetEnd
# will plot on the same figure the TDI signal for each phage. 
# if xSacle is True, x axis is scaled to relative position along the genome to allow comparrison of different length sequences. 
# if subset is true, each sequence is shortened to the positions defined in the final two fields of the input.
# other arguments (windowSize, stepSize, k) are passed to the TDI function. 
# plot titled with title. saves the resulting plot to saveName
# only first maxNum lines are plotted. the first 10 phage defined in the input file are plotted by default. 
# if saveFile is not None, saves a csv file of the data to the file specified. 
def compareTDI(nameFile, xScale, subset, windowSize, stepSize, k, title, saveName, maxNum, saveFile):
	# read information in nameFile
	names = []
	fnames = []
	subsets = []
	with open(nameFile, 'r') as nf:
		line = nf.readline().strip().split(',')
		while line != ['']:
			if subset and (len(line) < 4):
				print "Specified to subset but didn't include start and end position"
				print "Subsetting will be disabled"
				subset = False
			names.append(line[0])
			fnames.append(line[1])
			# parse subset information
			if subset:
				subsets.append([line[2],line[3]])
			line = nf.readline().strip().split(',')
	# compute for each defined phage
	data = []
	if subset:
		calcNum=0
		for fname, name, subset in zip(fnames, named, subsets):
			if calcNum<maxNum:
				print "Computing for " + name
				d = TDI(fname, windowSize, stepSize, k,  subset=subset)
				data.append(d)
			calcNum+=1
	else:
		calcNum=0
		for fname, name in zip(fnames, names):
			if calcNum<maxNum:
				print "Computing for " + name
				d = TDI(fname, windowSize, stepSize, k)
				data.append(d)
			calcNum+=1
	#if xScale, rescale axes and plot on genome position scale
	if xScale:
		for i in range(len(data)):
			scaledX = []
			for j in range(len(data[i][0])):
				scaledX.append(data[i][0][j] / float(data[i][0][len(data[i][0])-1]))
			data[i][0] = scaledX
	
	ax = plt.subplot(1,1,1)
	# Code from SO to make colors distinguishable
	if maxNum > len(names): numColors = len(names)
	else: numColors = numColors = maxNum

	cm = plt.get_cmap('gist_rainbow')
	cNorm  = colors.Normalize(vmin=0, vmax=numColors-1)
	scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
	ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(numColors)])

	i = 0
	for [xaxis, Zscores],name in zip(data, names):
		color = cm(1.*i/numColors)
		plt.plot(xaxis, Zscores, lw=1, label=name)
		i+=1

	handles, labels = ax.get_legend_handles_labels()
	plt.title(title)
	plt.xlabel('genomic position')
	plt.ylabel('TDI Z-score')
	plt.legend(prop={'size':5})
	plt.savefig(saveName, dpi=300)
	plt.clf()
	print "Plot \' " + title+ "\' saved to " + saveName 

	if saveFile != None:
		with open(saveFile, 'w') as sf:
			writer = csv.writer(sf)
			writer.writerow(['name']+data[0][0])
			for row, name  in zip(data, names):
				writer.writerow([name]+row[1])

	print "Done  :)"

# function call
compareTDI(nameFile, xScale, subset,windowSize,stepSize,k,title,saveName,maxNum,saveFile)