library(ape)
setwd("~/GitHub/tango")
metadata <- as.data.frame(read.table('data/phagesDB/sequenced_names_clusters.txt', header=T))
dev.4 <- as.matrix(read.table('data/kmer_counts/663_phage/all_dev_k4.tsv')) 
dev.4.clusters <- dev.4
rownames(dev.4.clusters) <- lapply(rownames(dev.4.clusters), function(x) strsplit(strsplit(x,'\\(')[[1]][2],'\\)')[[1]][1])

dev.4.clusters.log <- log2(dev.4.clusters+0.00001)

dist.4 <- dist(dev.4.clusters)
nj.4 <- nj(dist.4)
plot(nj.4)

dist.4.log <- dist(dev.4.clusters.log)
nj.4.log <- nj(dist.4.log)
plot(hclust(dist.4.log))

#trying pvclust package
library(pvclust)
result <- pvclust(t(dev.4.clusters), method.dist="cor", method.hclust="average", nboot=100)
plot(result)

dev.4.subset<- as.matrix(read.table('data/kmer_counts/Hatfull_60_subset/Hatful_60_dev_k4.tsv'))
res.subset <- pvclust(t(dev.4.subset))

seplot(res.subset, identify=TRUE)
# 15 29 are high SE