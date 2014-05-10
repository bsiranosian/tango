import matplotlib.pyplot as plt
from Bio import SeqIO
import math
import re
import numpy as np

## TOOLS FOR THE ANALYSIS OF TETRANUCLEOTIDE USAGE DEVAIATION IN GENOMES
## Designed primarily for phage genomes. May not be the fastest way, but it gets things done!

# reverse_complement(pattern): returns the reverse complement of DNA string pattern. Equivalent to an inversion. 
def reverse_complement(pattern):
	chars = list(pattern)
	complement = ''
	lookup = dict({("A","T"),("T","A"),("C","G"),("G","C"),("N","N")})
	for base in chars:
		complement += lookup[base]
	rev_complement = complement[::-1]
	return rev_complement 

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

#TUD(filename, k, kmerlist): computes the k-mer usage deviation for a genome in a fasta file defined in filename.
# also takes as input an integer k and a list of all k-mers over the nucleotide alphabet, generated by enumerateKmers 
# returns a dictionary of kmers as keys and normalized frequencies as values
# normalization done with a zero-order markov model (weighted frequencies)
def TUD(filename, k, kmerList, RC=False):
	#construct a dictionary of all kmers
	kmerDict = dict((key, 0 ) for key in kmerList)
	#parse fasta from filename
	for seq_record in SeqIO.parse(filename, "fasta"):
		sequence = seq_record.seq.tostring().upper()

	## NEW CODE ## 
	# append reverse complement if RC=True
	if RC:
		sequence += reverse_complement(sequence)

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

#similar to TUD(), but accepts a string instead of a filename
def TUDFromString(sequence, k, kmerList):
	kmerDict = dict((key, 0 ) for key in kmerList)
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

#similar to phageTUDCalc, except uses BioPython to parse through a single FASTA file with
#a large number of mycobacteria 
#modified to produce CSV files to reduce R processing errors
def mycobacteriaTUDCalc(filename, outfile, k):

	with open(outfile, 'w') as of:
		kmerList = enumerateKmers(4)
		wroteKmersAtTop = False

		for seq_record in SeqIO.parse(filename, "fasta"):
			sequence = seq_record.seq.tostring().upper()
			name = seq_record.description.split("| ")[1].replace(" ", "_").replace(',', '-')

			tud = TUDFromString(sequence, k, kmerList)

			#if it's the first time, write kmers at the top of the file
			#else, skip over
			if not wroteKmersAtTop:
				kmers = ''
				for kmer in tud.keys():
					kmers += '\t' + kmer
				of.write(kmers+'\n')
				wroteKmersAtTop = True

			print('working with ' + name)
			tud = TUDFromString(sequence, k, kmerList)
			#write each result to file
			toWrite = name
			for j in tud.values():
				toWrite += '\t' + str(j)
			of.write(toWrite+'\n')


#open a list of tab separated phage names (col1) and filenames (col2), one per line. 
#calculates TUD for each and writes results to outfile. 
def phageTUDCalc(phageFile, outfile, k,RC=False):
	#open and read in names and filenames
	with open(phageFile, 'r') as pf:
		phages = []
		line = pf.readline()
		while line != '':
			phages.append(line.strip().split('\t'))
			line=pf.readline()

	with open(outfile, 'w') as of:
		kmerList = enumerateKmers(k)
		#get one tud, write kmers at top
		oneTud = TUD(phages[0][1], k, kmerList,RC=RC)
		kmers = ''
		for kmer in oneTud.keys():
			kmers += '\t' + kmer
		of.write(kmers+'\n')
		#get TUD data for each phage
		for phage in phages:
			name = phage[0]
			print('working with ' + name)
			fname = phage[1]
			tud = TUD(fname, k, kmerList,RC=RC)
			#write each result to file
			toWrite = name
			for j in tud.values():
				toWrite += '\t' + str(j)
			of.write(toWrite+'\n')

#plots comparison line plot between the lines in dataFile (an output of phageTUDcalc)
# legend is the first row of dataFile. Title is displayed from title
# if plot is True, plots the resulting line graph. If false, returns the data structure of binned values.
# saves image to savename if plot is true. 
# plots frequency of log2 deviation. x scaled from -5 to 5
# if subset is defined, plot only the rows defined in the the list (0 indexed)
def plotTUD(dataFile, title, plot, saveName=None, subset=None):
	preNames = []
	preValues = []
	with open(dataFile, 'r') as df:
		line = df.readline().strip()
		#first line is the names of the kmers
		line = df.readline().strip()
		while line != '':
			preNames.append(line.split('\t')[0])
			preValues.append([math.log(float(v)+.00000001,2) for v in line.split('\t')[1:]])
			line = df.readline().strip()
	#if taking a subset, extract those
	values = []
	names =[]
	if subset != None:
		for i in subset:
			if i<len(preValues):
				values.append(preValues[i])
				names.append(preNames[i])
	else: 
		values = preValues
		names = preNames
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
			plt.plot(bins, binnedValues[i], lw=1, label=names[i])
		handles, labels = ax.get_legend_handles_labels()
		plt.title(title)
		plt.xlabel('log2(TUD)')
		plt.ylabel('frequency')
		plt.legend()
		plt.savefig(saveName, dpi=400)
		plt.clf()
	if not plot:
		return binnedValues

#TDI(): computes the k-mer difference index for a given genome. Takes in a filename, size of window to compute TDI in, 
# step size to move along the genome by, k, a list of kmers of length k, and a boolean plot. 
# if plot is true, plot the resulting figure using matplotlib. 
# always returns a list of window positions and Zscores
# if subset is defined, only return information about the sequence between the two positions. 
def TDI(filename, windowSize, stepSize, k, kmerList, subset=None, plot=False, RC=False):
	#construct a dictionary of all kmers
	kmerDict = dict((key, []) for key in kmerList)
	kmerDictGenome = dict((key, 0) for key in kmerList)
	#parse fasta from filename
	for seq_record in SeqIO.parse(filename, "fasta"):
		sequence = seq_record.seq.tostring().upper()
	# append reverse complement if RC=True
	if RC:
		sequence += reverse_complement(sequence)

	# pick out part defined in subset
	if subset != None:
		sequence = sequence[subset[0]:subset[1]]

	# slide along the genome in windows defined by windowSize with steps defined by stepSize
	start = 0
	end = windowSize
	windowIndex= []
	while end < len(sequence):
		if end > len(sequence):
			end = len(sequence)
		#print str(start) + " - " + str(end)
		window = sequence[start:end]
		#record list of windows
		windowIndex.append([start,end])
		for kmer in kmerDict.keys():
			oKmer, eKmer = TDI_kmerCount(kmer, window)
			# add observed over expected to dictionary
			if eKmer == 0:
				#print kmer
				#print oKmer
				kmerDict[kmer].append(1)
			else:
				kmerDict[kmer].append(oKmer/eKmer)
		#slide along the window by stepSize
		start += stepSize
		end += stepSize
	#compute kmers for whole genome.. not sure if we have to do this. 
	for kmer in kmerDictGenome.keys():
		oKmer, eKmer = TDI_kmerCount(kmer, sequence)
		# add observed over expected to dictionary
		kmerDictGenome[kmer]=oKmer/eKmer
	# compute tetranucleotide differences for each window. 
	differences = []
	for w in range(len(windowIndex)):
		# sum differences over all kmers
		sumD = 0
		for k, v in kmerDict.items():
			sumD += abs(kmerDictGenome[k] - v[w])
		differences.append(sumD)
	# compute Z-score   Z=(x-mu)/sigma
	Zscores = []
	mu = np.mean(differences)
	sigma = np.std(differences)
	for w in range(len(windowIndex)):
		Zscores.append((differences[w] - mu)/sigma)

	# PLOT SHIT
	xaxis = [x for x,y in windowIndex]
	if plot:
		plt.plot(xaxis, Zscores, lw=1)
		# ADD INFORMATION TO AXES, ET
		plt.show()

	return [xaxis,Zscores]

#helper function to count observed and expected kmers in a given window.
# takes in a kmer and window, both strings. Normalizes with markov chain method. 
def TDI_kmerCount(kmer, window):
	#observed count in window
	reg =  r'(?=('+kmer+'))'
	oKmer = float(len(re.findall(reg, window,)))
	#expected count in sequence done with markov chain analysis? Different than TUD calculation. 
	#ratio compared to the subwords. 
	reg1 = r'(?=('+kmer[1:]+'))'
	reg2 = r'(?=('+kmer[:len(kmer)-1]+'))'
	reg3 = r'(?=('+kmer[1:len(kmer)-1]+'))'
	#print "window: " + str(len(window))
	eKmer = (float(len(re.findall(reg1, window)))) * (float(len(re.findall(reg2, window,)))) / (float(len(re.findall(reg3, window))))
	#print (float(len(re.findall(reg1, window))))
	#print (float(len(re.findall(reg2, window)))) 
	#print (float(len(re.findall(reg3, window)))) 
	#print kmer + ": o: " +str(oKmer) + " e: "+str(eKmer)
	return (oKmer, eKmer)

#kmerCount(filename, k, kmerlist): computes the number of each kmer for a sequence. 
# also takes as input an integer k and a list of all k-mers over the nucleotide alphabet, generated by enumerateKmers 
# returns a dictionary of kmers as keys and counts as values
# if probability is true, returns the count / (n-k+1) as values
def kmerCount(sequence, k, kmerList, probability=False):
	#construct a dictionary of all kmers
	kmerDict = dict((key, 0 ) for key in kmerList)

	# slide a window of length k across the genome and add to dict
	for i in range(len(sequence)-k+1):
		if kmerDict.has_key(sequence[i:i+k]):
			kmerDict[sequence[i:i+k]] += 1
	# divide by number of kmers in the genome 
	if probability:
		for key in kmerDict.keys():
			kmerDict[key] = float(kmerDict[key])/(len(sequence)-k+1)
	return kmerDict

#calls kmerCount on a filename. see main function for details. 
# if RC is true, extends the sequcene by the reverse complement before counting
def doKmerCount(filename, k, probability=False, RC=False):
	kmerList = enumerateKmers(k)
	#parse fasta from filename
	for seq_record in SeqIO.parse(filename, "fasta"):
		sequence = seq_record.seq.tostring().upper()

	# append reverse complement if RC=True
	if RC:
		sequence += reverse_complement(sequence)
	return kmerCount(sequence, k, kmerList, probability)

#doKmerCountWindows: computes kemrCount at sliding windows across the genome. 
# returns (list of dictionaries one for each window, list of window start and end positions)
def doKmerCountWindows(filename, k, windowSize, stepSize, probability=False):
	kmerList = enumerateKmers(k)
	#parse fasta from filename
	for seq_record in SeqIO.parse(filename, "fasta"):
		sequence = seq_record.seq.tostring().upper()

	windows = int(math.ceil(len(sequence)/float(windowSize)))

	start = 0 
	end = start+windowSize
	toReturn = []
	while end < len(sequence):
		#append data
		toReturn.append((kmerCount(sequence[start:end], k, kmerList, probability),[start,end]))
		# get start and end positions
		start += stepSize
		end += stepSize
	#do one last calc for the last window 
	toReturn.append((kmerCount(sequence[start:end], k, kmerList, probability),[start,end]))
	return toReturn

#GCcontent: computes GC content in a sliding window across the genome. 
# returns (list of gc percentages, list of window starts, average GC content)
def GCcontent(sequence, windowSize, stepSize):
	start = 0
	end = start+windowSize
	starts = []
	GCpercents = []
	while end < len(sequence):
		sub = sequence[start:end]
		GCpercents.append((sub.count('G') + sub.count('C'))/float(len(sub)))
		starts.append(start)

		start+=stepSize
		end += stepSize
	#do one last calc for the last window
	GCpercents.append((sub.count('G') + sub.count('C'))/float(len(sub)))
	starts.append(start)

	return(GCpercents, starts, (sequence.count('G') + sequence.count('C'))/float(len(sequence)))

#doGCcontent: computes GC content for a sequence in a fasta file.
# see main function for details
# if plot=True, plots the resulting line graph
# if RC=True, passes the reverse complement of the sequence to the GC content function
def doGCcontent(filename, windowSize, stepSize, plot=False, RC=False):
	for seq_record in SeqIO.parse(filename, "fasta"):
		sequence = seq_record.seq.tostring().upper()

	if RC:
		sequence = reverse_complement(sequence)
	#do calc
	gc = GCcontent(sequence, windowSize, stepSize)
	if plot:
		plt.plot(gc[1],gc[0])
		plt.xlabel('genomic position')
		plt.ylabel('GC fraction')
		plt.show()
	return gc