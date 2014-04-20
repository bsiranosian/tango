#parameters 
import argparse

parser = argparse.ArgumentParser(description="This script produces a Tetranucleotide Difference Index (TDI) plot for multiple phage genomes. The only required inputs are the location of a configuration file, location to save the resulting figure and a title for the plot. The configuration file is a tab delimited file of phage info with 2 necessary fields and 2 optional fields: \n name\tfastaPath\tsubsetStart\tsubsetEnd\nEach phage to be compared should be separated by a new line. Additional arguments provide more control over the configuration of the plot and the underlying calculations.")
parser.add_argument("nameFile", action="store", metavar="nameFile", help="The file name of the tab delimited file defining phage information")
parser.add_argument("saveName", action="store", metavar="saveName", help="The file to save the resulting image to.")
parser.add_argument("title", action="store", metavar="title", help="The title for the resulting plot.")
parser.add_argument("--xScale", action="store", metavar="xScale", default="False", help="Set to True to scale the x axis of the plot to relative genome position.")
parser.add_argument("--subset", action="store", metavar="subset", default="False", help="set to True to plot only the regions defined in the input file. If True, each line must have integers in fields 3 and 4 that represent the genomic region of each phage to compare")
parser.add_argument("--windowSize", action="store", metavar="windowSize", default="5000", help="The size of the window to compute TDI within. Default of 5000bp is used in Pride et. al 2006")
parser.add_argument("--stepSize", action="store", metavar="stepSize", default="1000", help="How many bases to move the window along the genome at each iteration. Default of 1000bp is used in Pride et. al 2006")
parser.add_argument("--k", action="store", metavar="k", default="4", help="Can also use this to compute 2,3,5-mers, etc. UNTESTED!")
args=parser.parse_args()

nameFile=args.nameFile
saveName=args.saveName
title=args.title
if args.xScale == "True": xScale=True
else: xScale=False
if args.subset == "True": subset=True
else: subset=False
windowSize=int(args.windowSize)
stepSize=int(args.stepSize)
k=int(args.k)


from tetranucleotideAnalysis import *
#compareTDI(): plots multiple tetranucleotide difference index signals on the same plot. 
# takes as input a tab delimited file of phage info with 2 necessary fields and 2 optional fields:
# 	name 	fastaPath	subsetStart	subsetEnd
# will plot on the same figure the TDI signal for each phage. 
# if xSacle is True, x axis is scaled to relative position along the genome to allow comparrison of different length sequences. 
# if subset is true, each sequence is shortened to the positions defined in the final two fields of the input.
# other arguments (windowSize, stepSize, k) are passed to the TDI function. 
# plot titled with title. saves the resulting plot to saveName
def compareTDI(nameFile, xScale, subset, windowSize, stepSize, k, title, saveName):
	# read information in nameFile
	names = []
	fnames = []
	subsets = []
	with open(nameFile, 'r') as nf:
		line = nf.readline().strip().split('\t')
		while line != ['']:
			if subset and (len(line) < 4):
				print "Specified to subset but didn't include start and end position"
				print "Subsetting will be disabled"
				subset = False
			names.append(line[0])
			fnames.append(line[1])
			if subset:
				subsets.append([line[2],line[3]])
			line = nf.readline().strip().split('\t')
	# compute for each defined phage
	data = []
	kmerList = enumerateKmers(k)
	if subset:
		for fname, name, subset in zip(fnames, named, subsets):
			print "Computing for " + name
			d = TDI(fname, windowSize, stepSize, k, kmerList, subset=subset)
			data.append(d)
	else:
		for fname, name in zip(fnames, names):
			print "Computing for " + name
			d = TDI(fname, windowSize, stepSize, k, kmerList)
			data.append(d)

	#if xScale, rescale axes and plot on genome position scale
	if xScale:
		for i in range(len(data)):
			scaledX = []
			for j in range(len(data[i][0])):
				scaledX.append(data[i][0][j] / float(data[i][0][len(data[i][0])-1]))
			data[i][0] = scaledX
	ax = plt.subplot(1,1,1)
	for [xaxis, Zscores],name in zip(data, names):
		plt.plot(xaxis, Zscores, lw=1, label=name)
	handles, labels = ax.get_legend_handles_labels()
	plt.title(title)
	plt.xlabel('genomic position')
	plt.ylabel('TDI Z-score')
	plt.legend()
	plt.savefig(saveName, dpi=300)
	plt.clf()
	print "Plot \' " + title+ "\' saved to " + saveName 
	print "Done  :)"

# function call
compareTDI(nameFile, xScale, subset,windowSize,stepSize,k,title,saveName)