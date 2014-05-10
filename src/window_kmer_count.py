# compute kmer counts for a list of files, save each result to an individual matrix file
from tetranucleotideAnalysis import *
import numpy as np 
import matplotlib.pyplot as plt

# get phage and fastas
nameFile = 'C:/Users/Admin/Documents/GitHub/tango/data/TDI_individual_clusters/windows/sequenced_phage_map_M1.txt'
#nameFile = 'C:/Users/Admin/Documents/GitHub/tango/data/phagesDB/sequenced_phage_map_windows.txt'
with open(nameFile,'r') as nf:
	names=[]
	fileNames=[]
	line = nf.readline().strip().split('\t')
	while line != ['']:
		names.append(line[0])
		fileNames.append(line[1])
		line = nf.readline().strip().split('\t')

# SET WINDOW AND STEP
windowSize=5000
stepSize=2500

#save all phage windows to one file
allOutFile = 'C:/Users/Admin/Documents/GitHub/tango/data/window_kmer_count/all_freqs_2000_k4.tsv'
with open(allOutFile,'w') as of:
	kmerList = enumerateKmers(4)
	kmerString = ''
	for kmer in kmerList:
		kmerString += kmer + '\t'
	of.write(kmerString+'\n')
	for name,fileName in zip(names, fileNames):
		print name
		data=doKmerCountWindows(fileName, 4, windowSize, stepSize, probability=False)
		for dataDict, window in data:
			toWrite = name+'_'+str(window[0])+':'+str(window[1])+'\t'
			for kmer in kmerList:
				toWrite += '\t'+str(dataDict[kmer])
			of.write(toWrite+'\n')

#save results here 
outFolder = 'C:/Users/Admin/Documents/GitHub/tango/data/window_kmer_count/'
for name,fileName in zip(names, fileNames):
	with open(outFolder+name+'_'+str(windowSize)+'_'+str(stepSize)+ '_k4_prob' +'.txt','w') as of:
		data=doKmerCountWindows(fileName, 4, windowSize, stepSize, probability=True)
		# write header
		kmerList = enumerateKmers(4)
		kmerString = ''
		for kmer in kmerList:
			kmerString += kmer + '\t'
		of.write(kmerString+'\n')

		for dataDict, window in data:
			toWrite = str(window[0])+':'+str(window[1])+'\t'
			for kmer in kmerList:
				toWrite += '\t'+str(dataDict[kmer])
			of.write(toWrite+'\n')


#Do with di and trinucleotides
#save results here 
outFolder = 'C:/Users/Admin/Documents/GitHub/tango/data/window_kmer_count/'
for name,fileName in zip(names, fileNames):
	with open(outFolder+name+'_'+str(windowSize)+'_'+str(stepSize)+ '_k2' +'.txt','w') as of:
		data=doKmerCountWindows(fileName, 2, windowSize, stepSize, probability=True)
		# write header
		kmerList = enumerateKmers(2)
		kmerString = ''
		for kmer in kmerList:
			kmerString += kmer + '\t'
		of.write(kmerString+'\n')

		for dataDict, window in data:
			toWrite = str(window[0])+':'+str(window[1])+'\t'
			for kmer in kmerList:
				toWrite += '\t'+str(dataDict[kmer])
			of.write(toWrite+'\n')

outFolder = 'C:/Users/Admin/Documents/GitHub/tango/data/window_kmer_count/'
for name,fileName in zip(names, fileNames):
	with open(outFolder+name+'_'+str(windowSize)+'_'+str(stepSize)+ '_k3' +'.txt','w') as of:
		data=doKmerCountWindows(fileName, 3, windowSize, stepSize, probability=True)
		# write header
		kmerList = enumerateKmers(3)
		kmerString = ''
		for kmer in kmerList:
			kmerString += kmer + '\t'
		of.write(kmerString+'\n')

		for dataDict, window in data:
			toWrite = str(window[0])+':'+str(window[1])+'\t'
			for kmer in kmerList:
				toWrite += '\t'+str(dataDict[kmer])
			of.write(toWrite+'\n')
