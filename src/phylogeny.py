#working to calculate distance matrix -> phylogenetic tree
# distance between 2 phage is euclidean distance in 256 space
from tetranucleotideAnalysis import *

# calculates and saves the pairwise distance matrix between all phages
# reads from filenames, phage names stored in names, writes to outfile
def distance(filenames, names, outname,k):
	#build matrix of TUD numbers 
	TUDmatix = []
	kmerList = enumerateKmers(k)
	for f,n in zip(filenames,names):
		print 'Working with ' + n
		TUDmatix.append(TUD(f, k, kmerList).values())

	#fill in distanceMatrix
	print 'calculating distance matrix'
	distanceMatrix = [[0 for i in range(len(filenames))] for i in range(len(filenames))]
	for i in range(len(filenames)):
		print 'row ' + str(i) + ' of ' + str(len(filenames))
		for j in range(len(filenames)):
			dsum = 0
			for kmer in range(256):
				dsum += abs(TUDmatix[i][kmer] - TUDmatix[j][kmer])
			distanceMatrix[i][j] = (1.0/256.0) * dsum
	#write matrix to table
	fout = open(outname, 'w')
	#write names across top
	top = ''
	for name in names:
		top += ',' + name
	fout.write(top + '\n')
	for name, row in zip(names, distanceMatrix):
		toWrite = name
		for col in row:
			toWrite += ',' + str(col)
		fout.write(toWrite + '\n')

fdir = 'C:\\Users\\Admin\\Documents\\GitHub\\balboa\\phages\\fasta\\'
allPhageNames = ['Bethlehem','Bxb1','DD5','Jasper','KBG','Lockley','Solon','U2','Bxz2','Che12','D29','L5','Pukovnik','Chah','Orion','PG1','Qyrzula','Rosebush','Phaedrus','Pipefish','Cooper','Nigel','Bxz1','Cali','Catera','Rizal','ScottMcG','Spud','Myrna','Adjutor','Butterscotch','Gumball','PLot','PBI1','Troll4','244','Cjw1','Kostya','Porky','Boomer','Che8','Fruitloop','Llij','Pacc40','PMC','Ramsey','Tweety','Che9d','BPs','Halo','Konstantine','Predator','Barnyard','Brujita','Che9c','Corndog','Giles','Omega','TM4','Wildcat']
allPhageClusters = ['A1','A1','A1','A1','A1','A1','A1','A1','A2','A2','A2','A2','A2','B1','B1','B1','B2','B2','B3','B3','B4','B4','C1','C1','C1','C1','C1','C1','C2','D','D','D','D','D','D','E','E','E','E','F1','F1','F1','F1','F1','F1','F1','F1','F2','G','G','H1','H1','H2','I','I','Single','Single','Single','Single','Single']
allNames= [p + ' ('+c+')' for p,c in zip(allPhageNames, allPhageClusters)]
allFnames = [fdir+name+'.fasta' for name in allPhageNames]

distance(allFnames, allNames, 'allDistance.csv', 4)