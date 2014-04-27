from tetranucleotideAnalysis import *
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import numpy as np


#countMotif(sequence, motif, windowSize, stepSize): counts the occurances of a given motif in windows of windowSize, moving along sequence in intervals of stepSize
# returns a list of tuple: (position of start of region, occurances of the given motif for region) 
def countMotif(sequence, motif, windowSize, stepSize):
	start = 0
	end = start+windowSize
	counts =[]
	while end < len(sequence):
		window=sequence[start:end]
		#find all with regexp 
		reg =r'(?=('+motif+'))'
		count = len(re.findall(reg, window))
		counts.append((start, count))
		start += stepSize
		end += stepSize
	return counts


# code to plot the frequency of a given motif across the genome of multiple phage
def plotMotifGenome(nameFile, motif, windowSize, stepSize, title, saveName, maxNum, saveData=None, RC=False):
	# read information in nameFile
	names = []
	fnames = []
	with open(nameFile, 'r') as nf:
		line = nf.readline().strip().split('\t')
		while line != ['']:
			names.append(line[0])
			fnames.append(line[1])
			line = nf.readline().strip().split('\t')
	# do calculations for each sequence
	data=[]
	for fname, name in zip(fnames, names):
			#parse fasta from filename
		for seq_record in SeqIO.parse(fname, "fasta"):
			sequence = seq_record.seq.tostring().upper()

		## NEW CODE ## 
		# append reverse complement if RC=True
		if RC:
			sequence += reverse_complement(sequence)

		print "Computing for " + name
		d = countMotif(sequence, motif, windowSize, stepSize)
		data.append(d)

	# get from data into lists across genome
	

	ax = plt.subplot(1,1,1)
	# Code from SO to make colors distinguishable
	NUM_COLORS = maxNum
	cm = plt.get_cmap('gist_rainbow')
	cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
	scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
	ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])

	i = 0
	maxXaxis = []
	for phage,name in zip(data, names):
		xaxis = [start for (start, count) in phage]
		# find largest set of bins to save data with
		if len(xaxis) < len(maxXaxis):
			maxXaxis = xaxis
		counts  = [count for (start, count) in phage]
		if i < maxNum:
			color = cm(1.*i/NUM_COLORS)
			plt.plot(xaxis, counts, lw=1, label=name)
			i+=1

	handles, labels = ax.get_legend_handles_labels()
	plt.title(title)
	plt.xlabel('genomic position')
	plt.ylabel('count of motif: '+ motif)
	plt.legend(prop={'size':5})
	plt.savefig(saveName, dpi=300)
	plt.clf()
	print "Plot \' " + title+ "\' saved to " + saveName 
	
	if saveData != None:
		with open(saveData, 'w') as of:
			#save data to output file
			# header of motif and genomic positions
			header = motif
			for pos in maxXaxis:
				header += '\t'+str(pos)
			of.write(header+'\n')
			for phage,name in zip(data, names):
				counts  = [count for (start, count) in phage]
				toWrite = name
				for count in counts:
					toWrite += '\t'+str(count)
				of.write(toWrite+ '\n')
	print "Done  :)"


#plot this cool GATC motif for some clusters
# B3 is extremely elevated. look in comparison to a few others 
plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B3.txt' , 'GATC', 1000,900,'GATC frequency in B3 genomes. 1000/900', '../figures/GATC_motif_B3.png',10,saveData='../data/GATC_motif_B3.tsv')
plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B2.txt' , 'GATC', 1000,900,'GATC frequency in B2 genomes. 1000/900', '../figures/GATC_motif_B2.png',10,saveData='../data/GATC_motif_B2.tsv')
plotMotifGenome('../data/TDI_individual_clusters/windows/sequenced_phage_map_B1.txt' , 'GATC', 1000,900,'GATC frequency in B1 genomes. 1000/900', '../figures/GATC_motif_B1.png',10,saveData='../data/GATC_motif_B1.tsv')

