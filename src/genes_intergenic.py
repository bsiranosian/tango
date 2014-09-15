# Examine if tetranucleotide usage is the same in genes and intergenic regions
from Bio import Entrez
 
# ref.: http://wilke.openwetware.org/Parsing_Genbank_files_with_Biopython.html
 
# replace with your real email (optional):
Entrez.email = 'whatever@mail.com'
# accession id works, returns genbank format, looks in the 'nucleotide' database:
genes=Entrez.efetch(db='nucleotide',id='DQ398052.2',rettype='gb')
fasta=Entrez.efetch(db='nucleotide',id='DQ398052.2',rettype='fasta')
# parse out coordinates of genes
a= [b.strip() for b in genes.readlines()]
te = [t.split()[1] for t in a if t[0:3]=='CDS']
te2 = [t.split('(')[1].split(')')[0] if t[0]=='c' else t for t in te]
coordinates=[[int(t.split('..')[0]),int(t.split('..')[1])] for t in te2]

# parse out fasta sequence
seqLines=[f.strip() for f in fasta.readlines() if (f[0] != '>' and f[0]!='\n') ]
seq = ''.join(map(str, seqLines))

# get sequences from coordinates
geneSeqs = [seq[s[0]-1:s[1]] for s in coordinates]

# sequences in between genes - NOT VALID BECAUSE OF OVERLAPS
intergenicCoordinates = [[0, coordinates[0][0]-1]]
for c in range(len(coordinates)-1):
	intergenicCoordinates.append([coordinates[c][1]+1, coordinates[c+1][0]-1])
intergenicCoordinates.append([coordinates[len(coordinates)-1][1]+1, len(seq)])

intergenicSeqs =[seq[s[0]-1:s[1]] for s in intergenicCoordinates]

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

def kmerCount(sequence, k, probability=False):
	#construct a dictionary of all kmers
	kmerDict = dict((key, 0 ) for key in enumerateKmers(k))

	# slide a window of length k across the genome and add to dict
	for i in range(len(sequence)-k+1):
		if kmerDict.has_key(sequence[i:i+k]):
			kmerDict[sequence[i:i+k]] += 1
	# divide by number of kmers in the genome 
	if probability:
		for key in kmerDict.keys():
			kmerDict[key] = float(kmerDict[key])/(len(sequence)-k+1)
	return kmerDict

def zeroOrderExpected(sequence, k):
	#construct a dictionary of all kmers
	kmerDict = dict((key, 0 ) for key in enumerateKmers(k))
	nucleotides = kmerCount(sequence, 1, probability=True)

	# calc expeccted number for each kmer in the dict
	for kmer in kmerDict.keys():
		kA = kmer.count('A')
		kT = kmer.count('T')
		kC = kmer.count('C')
		kG = kmer.count('G')
		eKmer = (math.pow(nucleotides['A'],kA) * \
				 math.pow(nucleotides['T'],kT) * \
				 math.pow(nucleotides['C'],kC) * \
				 math.pow(nucleotides['G'],kG))* \
				(len(sequence)-k+1)
		kmerDict[kmer] = eKmer
	return kmerDict

# look at tud in each gene... is it consistent?
import matplotlib.pyplot as plt
# check length distribution
plt.hist([a[1]-a[0] for a in coordinates], bins=20)

# only take genes longer than 500bp
obs_4 = [kmerCount(a, 4) for a in geneSeqs if len(a) >499]
exp_4 = [zeroOrderExpected(a, 4) for a in geneSeqs if len(a) >499]
dev_4 = [[float(a[kmer])/float(b[kmer]) for kmer in a.keys()] for a,b in zip(obs_4, exp_4)]
# unfortunately lots of things are zero 
# try with trinucleotides
obs_3 = [kmerCount(a, 3) for a in geneSeqs if len(a) > 499]
exp_3 = [zeroOrderExpected(a, 3) for a in geneSeqs if len(a) > 499]
dev_3 = [[float(a[kmer])/float(b[kmer]) for kmer in a.keys()] for a,b in zip(obs_3, exp_3)]

with open('C:/Users/Admin/Documents/GitHub/tango/data/gene_analysis/wlidcat_gene_dev3.txt', 'w') as of:
	for d in dev_3:
		of.write(','.join([str(b) for b in d]) + '\n')

with open('C:/Users/Admin/Documents/GitHub/tango/data/gene_analysis/wlidcat_gene_dev4.txt', 'w') as of:
	for d in dev_4:
		of.write(','.join([str(b) for b in d]) + '\n')