library(ape)
setwd("~/GitHub/tango")
metadata <- as.data.frame(read.table('data/phagesDB/sequenced_names_clusters.txt', header=T))
dev.4 <- as.matrix(read.table('data/kmer_counts/663_phage/all_dev_k4.tsv')) 
dev.4.clusters <- dev.4
rownames(dev.4.clusters) <- lapply(rownames(dev.4.clusters), function(x) strsplit(strsplit(x,'\\(')[[1]][2],'\\)')[[1]][1])

dev.4.clusters.log <- log2(dev.4.clusters)

dist.4 <- dist(dev.4.clusters)
nj.4 <- nj(dist.4)
plot(nj.4)

dist.4.log <- dist(dev.4.clusters.log)
nj.4.log <- nj(dist.4.log)
plot(nj.4.log)
