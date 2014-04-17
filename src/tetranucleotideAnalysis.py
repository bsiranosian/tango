import matplotlib.pyplot as plt
from Bio import SeqIO
import math
import re

## TOOLS FOR THE ANALYSIS OF TETRANUCLEOTIDE USAGE DEVAIATION IN GENOMES
## Designed primarily for phage genomes. May not be the fastest way, but it gets things done!

#construct list of all kmers over alphabet ATCG
def enumerateKmers(k):
	alphabet = ['A','T','C','G']
	kmerList = ['']
	while k > 0:
		working = kmerList[:]
		kmerList = []
		for kmer in working:
			for letter in alphabet:
				kmerList.append(kmer+letter)
		k -= 1 
	return kmerList

def TUD(filename, k, kmerList):
	#construct a dictionary of all kmers
	kmerDict = dict((key, 0 ) for key in kmerList)
	#parse fasta from filename
	for seq_record in SeqIO.parse(filename, "fasta"):
		sequence = seq_record.seq.tostring().upper()

	#for normalization need number of each nucleotide
	eA = sequence.count('A')/float(len(sequence))
	eT = sequence.count('T')/float(len(sequence))
	eC = sequence.count('C')/float(len(sequence))
	eG = sequence.count('G')/float(len(sequence))
	#calculate TUD for each kmer
	for kmer in kmerDict.keys():
		#print 'calculating for ' + kmer
		#observed count in genome
		reg =  r'(?=('+kmer+'))'
		oKmer = float(len(re.findall(reg, sequence)))
		#print 'observed: ' + str(oKmer)
		#expected count in sequence
		#number of each letter in kmer 
		kA = kmer.count('A')
		kT = kmer.count('T')
		kC = kmer.count('C')
		kG = kmer.count('G')
		eKmer = (math.pow(eA,kA)*math.pow(eT,kT)*math.pow(eC,kC)*math.pow(eG,kG))*len(sequence)
		#frequency is observed/expected
		kmerDict[kmer] =oKmer/eKmer
	return kmerDict

#open a list of tab separated phage names (col1) and filenames (col2), one per line. 
#calculates TUD for each and writes results to outfile. 
def phageTUDCalc(phageFile, outfile, k):
	#open and read in names and filenames
	with open(phageFile, 'r') as pf:
		phages = []
		line = pf.readline()
		while line != '':
			phages.append(line.strip().split('\t'))
			line=pf.readline()

	with open(outfile, 'w') as of:
		kmerList = enumerateKmers(4)
		#get one tud, write kmers at top
		oneTud = TUD(phages[0][1], k, kmerList)
		kmers = ''
		for kmer in oneTud.keys():
			kmers += '\t' + kmer
		of.write(kmers+'\n')
		#get TUD data for each phage
		for phage in phages:
			name = phage[0]
			print 'working with ' + name
			fname = phage[1]
			tud = TUD(fname, k, kmerList)
			#write each result to file
			toWrite = name
			for j in tud.values():
				toWrite += '\t' + str(j)
			of.write(toWrite+'\n')

#plots comparison histograms between the fastas in filenames
#legend is titles with names provided in names
# if plot is True, plots the resulting line graph. If false, returns the data structure of binned values.
def plotHist(filenames, names, k, plot):
	values = []
	kmerList = enumerateKmers(k)
	for f, n in zip(filenames, names):
		print 'Working with ' + n
		tud = TUD(f, k, kmerList)
		toAppend = []
		for v in tud.values():
			if v >0:
				toAppend.append(math.log(v,2))
			else: toAppend.append(0)
		values.append(toAppend)

	#assign to bins
	bins = [i/4.0 for i in range(-20, 21)]
	binnedValues = [[0 for i in range(len(bins))] for j in range(len(values))]

	#calculate number in each bin
	for f in range(len(values)):
		for b in range(len(bins)-1):
			binnedValues[f][b] = len([a for a in values[f] if (a<bins[b+1]) and (a >= bins[b])])

	#plot each on the same graph
	if plot:
		ax = plt.subplot(1,1,1)
		for i in range(len(values)):
			plt.plot(bins, binnedValues[i], lw=2, label=names[i])
		handles, labels = ax.get_legend_handles_labels()
		plt.legend()
		plt.show()
	if not plot:
		return binnedValues
