setwd("~/GitHub/tango/data/distance_windows")
library(gplots)
tud_663 <- as.matrix(read.table('663_tud.csv',sep=','))
athena_2000_500 <- as.matrix(read.table('athena_2000_500_tud.csv',sep=','))
# cluster infor from phagesdb
phageData <- as.data.frame(read.table('../../kmer_analysis/examples/PhagesDB_Data.txt',skip=1))
#unique tetranucleotides becasue of RC
unique <- c('AAAA', 'AAAT', 'AAAC', 'AAAG', 'AATA', 'AATT', 'AATC', 'AATG', 'AACA', 'AACT', 'AACC', 'AACG', 'AAGA', 'AAGT', 'AAGC', 'AAGG', 'ATAA', 'ATAT', 'ATAC', 'ATAG', 'ATTA', 'ATTC', 'ATTG', 'ATCA', 'ATCT', 'ATCC', 'ATCG', 'ATGA', 'ATGT', 'ATGC', 'ATGG', 'ACAA', 'ACAC', 'ACAG', 'ACTA', 'ACTC', 'ACTG', 'ACCA', 'ACCT', 'ACCC', 'ACCG', 'ACGA', 'ACGT', 'ACGC', 'ACGG', 'AGAA', 'AGAC', 'AGAG', 'AGTA', 'AGTC', 'AGTG', 'AGCA', 'AGCT', 'AGCC', 'AGCG', 'AGGA', 'AGGC', 'AGGG', 'TAAA', 'TAAC', 'TAAG', 'TATA', 'TATC', 'TATG', 'TACA', 'TACC', 'TACG', 'TAGA', 'TAGC', 'TAGG', 'TTAA', 'TTAC', 'TTAG', 'TTTC', 'TTTG', 'TTCA', 'TTCC', 'TTCG', 'TTGA', 'TTGC', 'TTGG', 'TCAC', 'TCAG', 'TCTC', 'TCTG', 'TCCA', 'TCCC', 'TCCG', 'TCGA', 'TCGC', 'TCGG', 'TGAC', 'TGAG', 'TGTC', 'TGTG', 'TGCA', 'TGCC', 'TGCG', 'TGGC', 'TGGG', 'CAAC', 'CAAG', 'CATC', 'CATG', 'CACC', 'CACG', 'CAGC', 'CAGG', 'CTAC', 'CTAG', 'CTTC', 'CTCC', 'CTCG', 'CTGC', 'CTGG', 'CCAC', 'CCTC', 'CCCC', 'CCCG', 'CCGC', 'CCGG', 'CGAC', 'CGTC', 'CGCC', 'CGCG', 'CGGC', 'GAAC', 'GATC', 'GACC', 'GAGC', 'GTAC', 'GTCC', 'GTGC', 'GCCC', 'GCGC', 'GGCC')

distances <- apply(athena_2000_500, 1, function(x) dist(rbind(x[x!=0],tud_663[1,][x!=0])))
for (phage in rownames(tud_663[2:663,])) {
  distances <- rbind(distances,apply(athena_2000_500, 1, function(x) dist(rbind(x[x!=0],tud_663[phage,][x!=0]))))
}
rownames(distances) <- rownames(tud_663)

distances_w0 <- apply(athena_2000_500, 1, function(x) dist(rbind(x,tud_663[1,])))
for (phage in rownames(tud_663[2:663,])) {
  distances_w0 <- rbind(distances_w0,apply(athena_2000_500, 1, function(x) dist(rbind(x,tud_663[phage,]))))
}
rownames(distances_w0) <- rownames(tud_663)


#minimum distance
which(distances==min(distances),arr.ind = T)
# Kamiyu(B3)
plot(distances["Kamiyu(B3)",], col='red', type='o')
points(distances["Athena(B3)",], col='blue', type='o')

pdf('Athena_2000_500_distance.pdf', width=20, height=20)
heatmap.2(distances, Colv = NA, trace ="none")
dev.off()

pdf('Athena_2000_500_distance_w0.pdf', width=20, height=20)
heatmap.2(distances_w0, Colv = NA, trace ="none")
dev.off()

# including zeros in calculation
for (name in rownames(tud_663)){
  print(name)
  windows <- as.matrix(read.table(paste('2000_500_data/',name,'.csv',sep = ''),sep=','))
  
  distances_w0 <- apply(windows, 1, function(x) dist(rbind(x,tud_663[1,])))
  for (phage in rownames(tud_663[2:663,])) {
    distances_w0 <- rbind(distances_w0,apply(windows, 1, function(x) dist(rbind(x,tud_663[phage,]))))
  }
  rownames(distances_w0) <- rownames(tud_663)
  
  write.csv(distances_w0,paste('2000_500_data/',name,'_distances_with0.csv',sep = ''))

  pdf(paste('2000_500_heatmaps_with0/',name,'_with0.pdf',sep = ''), width=20, height=20, pointsize = 6)
  heatmap.2(distances_w0, Colv = NA, trace ="none", col=redblue(256))
  dev.off()
  
}

# without including zeros in calculation 
for (name in rownames(tud_663)){
  print(name)
  windows <- as.matrix(read.table(paste('2000_500_data/',name,'.csv',sep = ''),sep=','))
  
  distances <- apply(windows, 1, function(x) dist(rbind(x[x!=0],tud_663[1,][x!=0])))
  for (phage in rownames(tud_663[2:663,])) {
    distances <- rbind(distances,apply(windows, 1, function(x) dist(rbind(x[x!=0],tud_663[phage,][x!=0]))))
  }
  rownames(distances) <- rownames(tud_663)
  
  write.csv(distances,paste('2000_500_data/',name,'_distances.csv',sep = ''))
  
  pdf(paste('2000_500_heatmaps/',name,'.pdf',sep = ''), width=20, height=20, pointsize = 6)
  heatmap.2(distances, Colv = NA, trace ="none", col=redblue(256))
  dev.off()
  
}

# without including zeros in calculation - LOG 
for (name in rownames(tud_663)){
  print(name)
  windows <- as.matrix(read.table(paste('2000_500_data/',name,'.csv',sep = ''),sep=','))
  
  distances_log <- apply(windows, 1, function(x) dist(rbind(log(x[x>0 & tud_663[1,]>0]),log(tud_663[1,][x>0 & tud_663[1,]>0]))))
  for (phage in rownames(tud_663[2:663,])) {
    distances_log <- rbind(distances_log,apply(windows, 1, function(x) dist(rbind(log(x[x>0 & tud_663[phage,]>0]),log(tud_663[phage,][x>0 & tud_663[phage,]>0])))))
  }
  rownames(distances_log) <- rownames(tud_663)
  
  write.csv(distances_log,paste('2000_500_data/',name,'_distances_log.csv',sep = ''))
  
  pdf(paste('2000_500_heatmaps_log/',name,'.pdf',sep = ''), width=20, height=20, pointsize = 6)
  heatmap.2(distances_log, Colv = NA, trace ="none", col=redblue(256))
  dev.off()
  
}

# possible J to A11 (and other A subcluster) HGT at the beginning of A genomes
# compare two phage to in windows to see 
# def compare2(many) takes in a list of phage names, computes distances between
# windows of each. Returns distance matrix and plots
compare_many <- function(many){
  compare <- c()
  for (phage in many){
    compare <- rbind(compare,as.matrix(read.table(paste('2000_500_data/',phage,'.csv',sep = ''),sep=',')) )
  }
  # remove duplicate tetras - working with RC coutned data
  unique <- c('AAAA', 'AAAT', 'AAAC', 'AAAG', 'AATA', 'AATT', 'AATC', 'AATG', 'AACA', 'AACT', 'AACC', 'AACG', 'AAGA', 'AAGT', 'AAGC', 'AAGG', 'ATAA', 'ATAT', 'ATAC', 'ATAG', 'ATTA', 'ATTC', 'ATTG', 'ATCA', 'ATCT', 'ATCC', 'ATCG', 'ATGA', 'ATGT', 'ATGC', 'ATGG', 'ACAA', 'ACAC', 'ACAG', 'ACTA', 'ACTC', 'ACTG', 'ACCA', 'ACCT', 'ACCC', 'ACCG', 'ACGA', 'ACGT', 'ACGC', 'ACGG', 'AGAA', 'AGAC', 'AGAG', 'AGTA', 'AGTC', 'AGTG', 'AGCA', 'AGCT', 'AGCC', 'AGCG', 'AGGA', 'AGGC', 'AGGG', 'TAAA', 'TAAC', 'TAAG', 'TATA', 'TATC', 'TATG', 'TACA', 'TACC', 'TACG', 'TAGA', 'TAGC', 'TAGG', 'TTAA', 'TTAC', 'TTAG', 'TTTC', 'TTTG', 'TTCA', 'TTCC', 'TTCG', 'TTGA', 'TTGC', 'TTGG', 'TCAC', 'TCAG', 'TCTC', 'TCTG', 'TCCA', 'TCCC', 'TCCG', 'TCGA', 'TCGC', 'TCGG', 'TGAC', 'TGAG', 'TGTC', 'TGTG', 'TGCA', 'TGCC', 'TGCG', 'TGGC', 'TGGG', 'CAAC', 'CAAG', 'CATC', 'CATG', 'CACC', 'CACG', 'CAGC', 'CAGG', 'CTAC', 'CTAG', 'CTTC', 'CTCC', 'CTCG', 'CTGC', 'CTGG', 'CCAC', 'CCTC', 'CCCC', 'CCCG', 'CCGC', 'CCGG', 'CGAC', 'CGTC', 'CGCC', 'CGCG', 'CGGC', 'GAAC', 'GATC', 'GACC', 'GAGC', 'GTAC', 'GTCC', 'GTGC', 'GCCC', 'GCGC', 'GGCC')
  compare <- compare[,unique]
  compare[compare==0] <- NA
  d<-dist(compare) 
  dm <- as.matrix(d)
  #remove bins that overlap
  delta <- row(dm)-col(dm)
  dm[delta<4 & delta>-4] <- NA
  heatmap.2(dm, Colv = NA, Rowv=NA, trace ="none", col=redblue(256))
}


# OBSERVATIONS
# C1 - middle has really uncorrelated piece
# E - end has self-similar piece different than rest of genome
# is it similar to a singleton?

sings <- rownames(tud_663)[grep(pattern = '(Singleton)', rownames(tud_663))]
compare_list <- c('Tuco(E)', sings)
compare_many(compare_list)

#plot for a phage from each cluster
representative <- c()
for (name in names(table(phageData[,2]))){
  search <- paste('(',paste(name,')',sep=''),sep='')
  print(search)
  print(grep(search, rownames(tud_663),fixed=T, value = T)[1])
  representative <- c(representative, grep(search, rownames(tud_663),fixed=T, value = T)[1])
  representative <-representative[1:48]
}

for (name in representative){
  png(paste('2000_500_within_heatmaps/',name,'.png',sep = ''), width=2000, height=2000, pointsize = 10)
  compare_many(name)
  dev.off()
}

# H1 has a section in the middle much different than the rest of the genome. Is this a complexity issue
# repetitions at 35800
damien <- as.matrix(read.table(paste('2000_500_data/','Damien(H1)','.csv',sep = ''),sep=','))
damien <- damien[,unique]
num_zero <- apply(damien, 1, function(x) sum(x==0))

# F3 difference - 35-37k
