# doing Naive Bayes on tetranucleotide probabilities with multinomial distributions. 

#read in 4-mer probabilities and frequencies for all phage
probs <- as.data.frame(read.table('~/GitHub/tango/data/with_reverse_complement/FCGR_all_probability.tsv'))
freqs <- as.data.frame(read.table('~/GitHub/tango/data/with_reverse_complement/FCGR_all_frequency.tsv'))

#have to compute with integers (freqs)
#compute with JoeDirt against all others
compare <- freqs["JoeDirt(L1)",]
compared <- apply(probs, 1, function(x) dmultinom(compare, size=sum(compare),prob=x,log=T))
compared <- sort(compared,decreasing=T)
head(compared)
# Other L1s, then L2s, then L3 are the most likely. makes sense. 

#compare the last section of JoeDirt agains all others 
jd <- as.matrix(read.table('~/GitHub/tango/data/window_kmer_count//JoeDirt(L1)_2000_2000_k4_freq.txt'))
compare_start = apply(jd[1:32,], 2, sum)
compare_end = apply(jd[36:38,], 2, sum)
#compare_end=jd[38,]
compared_start = apply(probs, 1, function(x) dmultinom(compare_start, size=sum(compare_start),prob=x,log=T))
compared_end = apply(probs, 1, function(x) dmultinom(compare_end, size=sum(compare_end),prob=x,log=T))
compared_start <- sort(compared_start,decreasing=T)
compared_end <- sort(compared_end,decreasing=T)
head(compared_start,n=15)
head(compared_end,n=30)

#compare last of jd against all windows
probs <- as.data.frame(read.table('~/GitHub/tango/data/window_kmer_count/all_prob_5000_k4.tsv'))
freqs <- as.data.frame(read.table('~/GitHub/tango/data/window_kmer_count//all_freq_5000_k4.tsv'))
compare_end = freqs["JoeDirt(L1)_70000:75000",]
compared_end = apply(probs, 1, function(x) dmultinom(compare_end, size=sum(compare_end),prob=x,log=T))
compared_end <- sort(compared_end,decreasing=T)
head(compared_end,n=100)
