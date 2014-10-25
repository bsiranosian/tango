# Want to see if distance can be used to detect HGT better
# compute distance in a sliding window compared to all other phage

from kmerAnalysis import *
	

# readFastaMap(fastaMap): makes reading fastaMap easier
def readFastaMap(fastaMap):
	# read information in fastaMap
	names = []
	fnames = []
	subsets = []
	with open(fastaMap, 'r') as nf:
		line = nf.readline().strip().split(',')
		while line != ['']:
			names.append(line[0])
			fnames.append(line[1])
			# parse subset information
			line = nf.readline().strip().split(',')
	return [names,fnames]

# doTUD for all in fastaMap
def doTUD(fastaMap, k, RC):
	read = readFastaMap(fastaMap)
	names = read[0]
	fnames = read[1]
	
	# compute deviation for each defined phage
	devDict = dict()
	for fname, name in zip(fnames, names):
		obs = doKmerCount(fname,k, RC=RC)
		exp = doZeroOrderExpected(fname, k, RC=RC) 
		devDict[name] = {kmer: float(obs[kmer])/exp[kmer] for kmer in obs.keys()}
	return devDict

# doTUDsingleWindows(fname, k ,RC, windowSize, stepSize): computes TUD for the fasta file
# returns (list of dictionaries one for each window, list of window start and end positions)
def doTUDsingleWindows(fname, k, RC, windowSize, stepSize):
	obs = doKmerCountWindows(fname, k, windowSize, stepSize)
	exp = doZeroOrderExpectedWindows(fname, k ,windowSize, stepSize)

	#calculate obs over expect
	windowDevList=[]
	for window1,window2 in zip(obs, exp):
		obsDict=window1[0]
		pos1=window1[1]
		expDict=window2[0]
		pos2=window2[1]

		assert pos1==pos2

		windowDevList.append([{kmer : float(obsDict[kmer])/expDict[kmer] for kmer in obsDict.keys()}, pos1])

	return windowDevList

# commands to do analysis
os.chdir("C:/Users/Admin/Documents/GitHub/tango/kmer_analysis/")
mapfile='../data/fasta_maps/windows/sequenced_phage_map_663.txt'
fmap = readFastaMap(mapfile)
tud_663  = doTUD(mapfile, 4, True)

# do in windows for Athena(B3)
athena_2000_500=doTUDsingleWindows(fmap[1][fmap[0].index('Athena(B3)')], 4, True, 2000,500)

# export this shit to R
# save TUD 
with open('../data/distance_windows/663_tud.csv','w') as of:
	# write header of file - names of tetras
	of.write(','.join(tud_663['Stinger(B4)'].keys())+'\n')
	for phage in tud_663.keys():
		of.write(str(phage) + ',' + ','.join([str(tud_663[phage][key]) for key in tud_663[phage].keys()]) + '\n')

#save athena windows
with open('../data/distance_windows/athena_2000_500_tud.csv','w') as of:
	# write header of file - names of tetras
	of.write(','.join(athena_2000_500[0][0].keys())+'\n')
	for windowWrap in athena_2000_500:
		of.write(str(windowWrap[1][0]) + ':' + str(windowWrap[1][1]) + ',' + ','.join([str(windowWrap[0][key]) for key in windowWrap[0].keys()]) + '\n')

#do for all phage
for name, fname in zip(fmap[0], fmap[1]):
	data = doTUDsingleWindows(fname, 4, True, 2000,500)
	#save athena windows
	with open('../data/distance_windows/2000_500_data/' + name + '.csv','w') as of:
		# write header of file - names of tetras
		of.write(','.join(data[0][0].keys())+'\n')
		for windowWrap in data:
			of.write(str(windowWrap[1][0]) + ':' + str(windowWrap[1][1]) + ',' + ','.join([str(windowWrap[0][key]) for key in windowWrap[0].keys()]) + '\n')
