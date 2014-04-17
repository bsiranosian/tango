#generate figures from paper
from tetranucleotideAnalysis import *
fdir = 'C:\\Users\\Admin\\Documents\\GitHub\\balboa\\phages\\fasta\\'

# plot of average across all phages - using 60 from Hatful 2010 paper
def plot1():
	allNames = ['Bethlehem','Bxb1','DD5','Jasper','KBG','Lockley','Solon','U2','Bxz2','Che12','D29','L5','Pukovnik','Chah','Orion','PG1','Qyrzula','Rosebush','Phaedrus','Pipefish','Cooper','Nigel','Bxz1','Cali','Catera','Rizal','ScottMcG','Spud','Myrna','Adjutor','Butterscotch','Gumball','PLot','PBI1','Troll4','244','Cjw1','Kostya','Porky','Boomer','Che8','Fruitloop','Llij','Pacc40','PMC','Ramsey','Tweety','Che9d','BPs','Halo','Konstantine','Predator','Barnyard','Brujita','Che9c','Corndog','Giles','Omega','TM4','Wildcat']
	allFnames = [fdir+name+'.fasta' for name in allNames]
	allValues = plotHist(allFnames, allNames, 4, False)
	#save allValues to a file
	fout = open('allValuesOut.txt', 'w')
	for v, n in zip(allValues, allNames):
		toWrite = ''
		for val in v:
			toWrite = toWrite + ',' + str(val)
		fout.write(n+toWrite+'\n')
	fout.close()
	#calculate average across all phages
	bins = [i/4.0 for i in range(-20, 21)]
	phageMean = []
	for i in range(len(bins)):
		total = 0
		for j in range(len(allNames)):
			total += allValues[j][i]
		phageMean.append(total/len(allNames))

	#plot the mean number directly
	ax = plt.subplot(1,1,1)
	plt.plot(bins, phageMean, lw=2, label='mean across 60 phages')
	handles, labels = ax.get_legend_handles_labels()
	plt.legend()
	plt.show()	

# A representative phage from each cluster on the same plot
def plot2():
	clusterNames = ['Bethlehem','Chah','Bxz1','Adjutor','244','Boomer','BPs','Konstantine','Brujita','Corndog']
	clusters = ['A','B','C','D','E','F','G','H','I','Sin']
	clusterFnames = [fdir+name+'.fasta' for name in clusterNames]
	clusterLegend = [n + ' ('+c+')' for n, c in zip(clusterNames,clusters)]
	plotHist(clusterFnames, clusterLegend, 4, True)


# Similarity within a given cluster (A1)
def plot3():
	a1Names = ['Bethlehem','Bxb1','DD5','Jasper','KBG','Lockley','Solon','U2']
	cluster = 'A1'
	legend = [n + ' ('+cluster+')' for n in a1Names]
	a1Fnames = [fdir+name+'.fasta' for name in a1Names]
	plotHist(a1Fnames, legend, 4, True)
# Similarity within a given cluster (sin)
def plot4():
	sinNames = ['Corndog','Giles','Omega','Wildcat']
	cluster = 'sin'
	legend = [n + ' ('+cluster+')' for n in sinNames]
	sinFnames = [fdir+name+'.fasta' for name in sinNames]
	plotHist(sinFnames, legend, 4, True)

plot4()





